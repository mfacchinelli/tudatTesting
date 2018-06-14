/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Mooij, E. "Orbit-State Model Selection for Solar-Sailing Mission Optimization."
 *          AIAA/AAS Astrodynamics Specialist Conference. 2012.
 */

#include <ctime>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Mathematics/Statistics/basicStatistics.h"

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "propagatorsComparison.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
    if ( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if ( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

//! Function to compute the next step size [Mooij, E.]
std::pair< double, bool > computeNewStepSize(
        const double stepSize,
        const std::pair< double, double >& orders,
        const double safetyFactorForNextStepSize,
        const std::pair< double, double >& minimumAndMaximumFactorsForNextStepSize,
        const Eigen::VectorXd& relativeErrorTolerance,
        const Eigen::VectorXd& absoluteErrorTolerance,
        const Eigen::VectorXd& lowerOrderEstimate,
        const Eigen::VectorXd& higherOrderEstimate )
{
    // Compute local truncation error
    const Eigen::VectorXd localTruncationError = lowerOrderEstimate - higherOrderEstimate;

    // Compute norm of error
    double currentErrorTolerance = 0.0;
    double localTruncationErrorNorm = 0;
    for ( unsigned int i = 0; i < lowerOrderEstimate.rows( ); i++ )
    {
        currentErrorTolerance = absoluteErrorTolerance[ i ] / relativeErrorTolerance[ i ];
        localTruncationErrorNorm += std::pow( localTruncationError[ i ] / (
                                                  std::max( std::fabs( higherOrderEstimate[ i ] ), currentErrorTolerance ) ), 2 );
        // incorrect, also need value from next time step (?)
    }
    localTruncationErrorNorm = std::sqrt( localTruncationErrorNorm );

    // Compute new step size
    double newStepSize;
    const double relativeTruncationError = localTruncationErrorNorm / relativeErrorTolerance.minCoeff( );
    const double referenceValue = safetyFactorForNextStepSize * std::pow( 1.0 / relativeTruncationError, 1.0 / orders.second );
    if ( localTruncationErrorNorm > relativeErrorTolerance.minCoeff( ) )
    {
        newStepSize = std::max( minimumAndMaximumFactorsForNextStepSize.first, referenceValue ) * stepSize;
    }
    else
    {
        if ( referenceValue < minimumAndMaximumFactorsForNextStepSize.second )
        {
            newStepSize = referenceValue * stepSize;
        }
        else
        {
            newStepSize = minimumAndMaximumFactorsForNextStepSize.second * stepSize;
        }
    }

    // Check if the current state can be accepted.
    const bool isIntegrationStepAccepted = relativeTruncationError <= 1.0;

    // Return the computed new step size and whether the current step taken is acceptable. If it
    // isn't, the step will be recomputed with the computed new step size.
    return std::make_pair( newStepSize, isIntegrationStepAccepted );
}

//! Execute propagation of orbit of Satellite around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::ephemerides;
    using namespace tudat::input_output;
    using namespace tudat::aerodynamics;
    using namespace tudat::unit_conversions;
    using namespace tudat::basic_astrodynamics;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            OUTERMOST LOOP                //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Select case
    //      0: Aerocapture
    //      1: Full aerobraking
    //      2: Interplanetary trajectory
    //      3: Circular orbit at LEO (Low Earth Orbit)
    //      4: Molniya orbit
    //      5: Low-thrust trajectory
    const std::vector< std::string > pathAdditionTestCase = { "aero", "aero_full", "inter", "circ", "moln", "low_thrust" };
    const bool singleTestCase = false;
    unsigned int initialTestCase = 2;

    // Start loop
    for ( unsigned int testCase = initialTestCase; testCase < 6; testCase++ )
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            DEFINE TEST CASES             //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Predefine variables
        double simulationDuration = 0.0;
        std::vector< double > constantTimeStep;
        std::string simulationCentralBody;
        int limitingSphericalHarminics;
        Eigen::Vector6d satelliteInitialStateInKeplerianElements;

        // Variables that can be optionally overwritten
        int referenceTolerances = -15;
        unsigned int numberOfSimulationCases = 1;
        std::vector< unsigned int > stateSizes = { 6, 6, 6, 6, 7, 7, 7, 6 };

        // Set tolerances and other state related parameters
        std::vector< int > relativeTolerances = { -13, -12, -11, -10, -9, -8, -7 };

        bool absoluteTolerancesSet = false;
        double velocityTolerance = 1e-6;
        double quaternionTolerance = 1e-15;
        Eigen::VectorXd cartesianTolerances;
        Eigen::VectorXd usmTolerances;

        // Thrust values for 6th case
        std::vector< double > thrustMagnitudes;

        // Select conditions based on case
        switch ( testCase )
        {
        case 0: // Aerocapture
        {
            // Set simulation time settings
            simulationDuration = 0.625 * physical_constants::JULIAN_DAY;
            constantTimeStep = { 1.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0 };

            // Set simulation central body
            simulationCentralBody = "Mars";
            limitingSphericalHarminics = 21;

            // Initial conditions
            satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = -17305.0e3;
            satelliteInitialStateInKeplerianElements( eccentricityIndex ) = 1.2;
            satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
            satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
            satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 23.4 );
            satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 230.0 );
            break;
        }
        case 1: // Full aerobraking
        {
            // Set simulation time settings.
            simulationDuration = 145.0 * physical_constants::JULIAN_DAY;
            constantTimeStep = { 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0 };

            // Set simulation central body
            simulationCentralBody = "Mars";
            limitingSphericalHarminics = 21;

            // Initial conditions
            satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 27228500.0;
            satelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.869218;
            satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
            satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
            satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 23.4 );
            satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );
            break;
        }
        case 2: // Interplanetary trajectory
        {
            // Set simulation time settings.
            simulationDuration = 260.0 * physical_constants::JULIAN_DAY;
            constantTimeStep = { 25.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1200.0 };

            // Set simulation central body
            simulationCentralBody = "Sun";

            // Initial conditions
            satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 188768611.5e3;
            satelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.207505835788806;
            satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 2.5 );
            satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 328.6 );
            satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 147.8 );
            satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
            break;
        }
        case 3: // Vittaldev circular orbit at LEO
        {
            // Set simulation time settings.
            simulationDuration = 10.0 * physical_constants::JULIAN_DAY;
            constantTimeStep = { 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0 };

            // Set simulation central body
            simulationCentralBody = "Earth";
            limitingSphericalHarminics = 2;

            // Initial conditions
            satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6936.0e3;
            satelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.0;
            satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 28.5 );
            satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 194.8 );
            satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 272.3 );
            satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
            break;
        }
        case 4: // Vittaldev Molniya orbit
        {
            // Set simulation time settings.
            simulationDuration = 25.0 * physical_constants::JULIAN_DAY;
            constantTimeStep = { 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0 };

            // Set simulation central body
            simulationCentralBody = "Earth";
            limitingSphericalHarminics = 2;

            // Initial conditions
            satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 26559.0e3;
            satelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.70;
            satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 63.2 );
            satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 206.3 );
            satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 281.6 );
            satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
            break;
        }
        case 5: // Vittaldev low-thrust orbit
        {
            // Set simulation time settings.
            simulationDuration = 10.0 * physical_constants::JULIAN_DAY;
            constantTimeStep = { 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0 };

            // Set custom tolerances to include mass
            stateSizes = { 7, 7, 7, 7, 8, 8, 8, 7 };
            absoluteTolerancesSet = true;

            cartesianTolerances = Eigen::VectorXd::Zero( 7 );
            cartesianTolerances.segment( 0, 6 ) = Eigen::VectorXd::Constant( 6, velocityTolerance );
            cartesianTolerances[ 6 ] = 1;

            usmTolerances = Eigen::VectorXd::Zero( 8 );
            usmTolerances.segment( 0, 3 ) = Eigen::VectorXd::Constant( 3, velocityTolerance );
            usmTolerances.segment( 3, 4 ) = Eigen::VectorXd::Constant( 4, quaternionTolerance );
            usmTolerances[ 7 ] = 1;

            // Set simulation central body
            simulationCentralBody = "Earth";
            limitingSphericalHarminics = 21;

            // Set thurst settings
            thrustMagnitudes = { 5e-3, 5e-2, 5e-1, 5e0, 5e1 };
            numberOfSimulationCases = thrustMagnitudes.size( );

            // Initial conditions
            satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6378.1363e3 + 1000e3;
            satelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.0;
            satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 0.0 );
            satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 0.0 );
            satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 0.0 );
            satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
            break;
        }
        default:
        {
            throw std::runtime_error( "Mode not recognized." );
        }
        }

        // Fill in standard tolerances
        if ( !absoluteTolerancesSet )
        {
            cartesianTolerances = Eigen::VectorXd::Zero( 6 );
            cartesianTolerances = Eigen::VectorXd::Constant( 6, velocityTolerance );
            usmTolerances = Eigen::VectorXd::Zero( 7 );
            usmTolerances.segment( 0, 3 ) = Eigen::VectorXd::Constant( 3, velocityTolerance );
            usmTolerances.segment( 3, 4 ) = Eigen::VectorXd::Constant( 4, quaternionTolerance );
        }

        // Give error is folder is not defined
        if ( testCase > pathAdditionTestCase.size( ) )
        {
            throw std::runtime_error( "Output folder is undefined for this test case." );
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Load Spice kernels
        spice_interface::loadStandardSpiceKernels( );

        // Set simulation time settings
        const double simulationStartEpoch = 7.0 * physical_constants::JULIAN_YEAR + 30.0 * 6.0 * physical_constants::JULIAN_DAY;
        const double simulationEndEpoch = simulationDuration + simulationStartEpoch;

        // Define body settings for simulation
        std::vector< std::string > bodiesToCreate;
        switch ( testCase )
        {
        case 0:
        case 1:
        {
            bodiesToCreate.push_back( simulationCentralBody );
            bodiesToCreate.push_back( "Sun" );
            bodiesToCreate.push_back( "Earth" );
            break;
        }
        case 2:
        {
            bodiesToCreate.push_back( simulationCentralBody );
            bodiesToCreate.push_back( "Mercury" );
            bodiesToCreate.push_back( "Mars" );
            bodiesToCreate.push_back( "Earth" );
            bodiesToCreate.push_back( "Jupiter" );
            bodiesToCreate.push_back( "Saturn" );
            bodiesToCreate.push_back( "Neptune" );
            bodiesToCreate.push_back( "Uranus" );
            break;
        }
        case 3:
        case 4:
        case 5:
        {
            bodiesToCreate.push_back( simulationCentralBody );
            bodiesToCreate.push_back( "Sun" );
            bodiesToCreate.push_back( "Moon" );
            break;
        }
        }

        // Create body objects
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 1.0e3, simulationEndEpoch + 1.0e3 );
        for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
            bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        }
        switch ( testCase )
        {
        case 0:
        case 1:
        {
            // Tabulated atmosphere settings
            std::map< int, std::string > tabulatedAtmosphereFiles;
            tabulatedAtmosphereFiles[ 0 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/density.dat";
            tabulatedAtmosphereFiles[ 1 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/pressure.dat";
            tabulatedAtmosphereFiles[ 2 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/temperature.dat";
            tabulatedAtmosphereFiles[ 3 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/gasConstant.dat";
            tabulatedAtmosphereFiles[ 4 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";
            std::vector< AtmosphereDependentVariables > atmosphereDependentVariables = {
                density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere,
                gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere };
            std::vector< AtmosphereIndependentVariables > atmosphereIndependentVariables = {
                longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };
            std::vector< interpolators::BoundaryInterpolationType > boundaryConditions = {
                interpolators::use_boundary_value, interpolators::use_boundary_value, interpolators::use_default_value };

            // Define default extrapolation values
            std::vector< std::vector< std::pair< double, double > > > extrapolationValues =
                    std::vector< std::vector< std::pair< double, double > > >(
                        5, std::vector< std::pair< double, double > >( 3, std::make_pair( 0.0, 0.0 ) ) );
            extrapolationValues.at( 0 ).at( 2 ) = { 7.62e-5, 0.0 };
            extrapolationValues.at( 1 ).at( 2 ) = { 2.35, 0.0 };
            extrapolationValues.at( 2 ).at( 2 ) = { 1.61e2, 186.813 };
            extrapolationValues.at( 3 ).at( 2 ) = { 190.7, 8183.0 };
            extrapolationValues.at( 4 ).at( 2 ) = { 1.377, 1.667 };

            // Generate gravitational and atmospheric body settings
            bodySettings[ simulationCentralBody ]->gravityFieldSettings = boost::make_shared<
                    FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
            bodySettings[ simulationCentralBody ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                        tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables,
                        boundaryConditions, extrapolationValues );
            break;
        }
        case 2:
        {
            bodySettings[ simulationCentralBody ]->gravityFieldSettings = boost::make_shared<
                    CentralGravityFieldSettings >( 1.32712440019e20 );
            break;
        }
        case 3:
        case 4:
        case 5:
        {
            bodySettings[ simulationCentralBody ]->gravityFieldSettings = boost::make_shared<
                    FromFileSphericalHarmonicsGravityFieldSettings >( ggm02s );
            bodySettings[ simulationCentralBody ]->atmosphereSettings = boost::make_shared< ExponentialAtmosphereSettings >(
                        7050.0, 240.0, 1.225, 2.87e2 );
            break;
        }
        }
        NamedBodyMap bodyMap = createBodies( bodySettings );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create spacecraft object
        bodyMap[ "Satellite" ] = boost::make_shared< Body >( );
        const double satelliteMass = 1000.0;
        bodyMap[ "Satellite" ]->setConstantBodyMass( satelliteMass );

        // Aerodynamic parameters
        const bool readAerodynamicCoefficientsFromFile = true;
        const double referenceAreaAerodynamic = 37.5;
        boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;
        if ( readAerodynamicCoefficientsFromFile )
        {
            // Aerodynamic coefficients from file
            std::map< int, std::string > aerodynamicCoefficientFiles;
            aerodynamicCoefficientFiles[ 0 ] = getTudatRootPath( ) + "External/MRODragCoefficients.txt";
            aerodynamicCoefficientFiles[ 2 ] = getTudatRootPath( ) + "External/MROLiftCoefficients.txt";

            // Create aerodynamic coefficient interface settings.
            aerodynamicCoefficientSettings = readTabulatedAerodynamicCoefficientsFromFiles(
                        aerodynamicCoefficientFiles, referenceAreaAerodynamic,
                        boost::assign::list_of( aerodynamics::angle_of_attack_dependent )( aerodynamics::altitude_dependent ),
                        true, true );
        }
        else
        {
            // Set constant aerodynamic drag coefficient
            const Eigen::Vector3d aerodynamicCoefficients = 2.2 * Eigen::Vector3d::UnitX( );
            aerodynamicCoefficientSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                        referenceAreaAerodynamic, aerodynamicCoefficients, true, true );
        }
        bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );

        // Create radiation pressure settings
        double referenceAreaRadiation = referenceAreaAerodynamic;
        double radiationPressureCoefficient = 1.25;
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back( simulationCentralBody );
        boost::shared_ptr< RadiationPressureInterfaceSettings > SatelliteRadiationPressureSettings =
                boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Create and set radiation pressure settings
        bodyMap[ "Satellite" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface(
                                                                   SatelliteRadiationPressureSettings, "Satellite", bodyMap ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Loop over simulation cases
        for ( unsigned int simulationCase = 0; simulationCase < numberOfSimulationCases; simulationCase++ )
        {
            // Progress
            std::cout << std::endl << "Simulation " << simulationCase + 1 << " out of " << numberOfSimulationCases << std::endl;

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;

            // Define propagation settings.
            std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
            switch ( testCase )
            {
            case 5:
            {
                // Define thrust settings
                double specificImpulse = 5000.0;
                boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
                        boost::make_shared< ThrustDirectionFromStateGuidanceSettings >( simulationCentralBody, true, false );
                boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
                        boost::make_shared< ConstantThrustEngineSettings >( thrustMagnitudes.at( simulationCase ), specificImpulse );

                // Define thrust acceleration settings
                accelerationsOfSatellite[ "Satellite" ].push_back(
                            boost::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );
                // no break, so add more accelerations
            }
            case 0:
            case 1:
            case 3:
            case 4:
            {
                // Switch between Kepler orbit or perturbed environment
                bool keplerOrbit = false;
                if ( keplerOrbit )
                {
                    // Only central gravity
                    accelerationsOfSatellite[ simulationCentralBody ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
                }
                else
                {
                    // Define spherical harmonics, third bodies, solar radiation and aerodynamic forces
                    accelerationsOfSatellite[ simulationCentralBody ].push_back(
                                boost::make_shared< SphericalHarmonicAccelerationSettings >( limitingSphericalHarminics,
                                                                                             limitingSphericalHarminics ) );
                    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
                    {
                        if ( bodiesToCreate.at( i ) != simulationCentralBody )
                        {
                            accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >(
                                                                                              central_gravity ) );
                        }
                    }
                    accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
                    accelerationsOfSatellite[ simulationCentralBody ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
                }
                break;
            }
            case 2:
            {
                // Define central gravity for central and third bodies and solar radiation
                for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
                {
                    accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
                }
                accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
                break;
            }
            }

            // Add acceleration information
            accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
            bodiesToPropagate.push_back( "Satellite" );
            centralBodies.push_back( simulationCentralBody );

            AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             CREATE INITIAL CONDITIONS              ////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Convert to Cartesian elements
            double mainGravitationalParameter = bodyMap.at( simulationCentralBody )->getGravityFieldModel( )->getGravitationalParameter( );
            const Eigen::Vector6d SatelliteInitialState = convertKeplerianToCartesianElements(
                        satelliteInitialStateInKeplerianElements, mainGravitationalParameter );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             LOOP OVER PROPAGATORS                  ////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Define absolute tolerances
            std::vector< Eigen::VectorXd > absoluteTolerances( 8 );
            absoluteTolerances.at( 0 ) = cartesianTolerances;
            absoluteTolerances.at( 4 ) = usmTolerances;
            absoluteTolerances.at( 5 ) = usmTolerances;
            absoluteTolerances.at( 6 ) = usmTolerances;
            absoluteTolerances.at( 7 ) = cartesianTolerances;

            // Loop over propagators
            int integratorLimit, valueLimit;
            double relativeTolerance;
            std::vector< string > nameAdditionPropagator = { "_cowell", "_encke", "_kepl", "_equi", "_usm7", "_usm6", "_usmem", "_ref" };
            std::vector< string > nameAdditionIntegrator = { "_var", "_const" };
            for ( int propagatorType = 0; propagatorType < 8; propagatorType++ )
            {
                // Skip Encke and Gauss propagators
                if ( !( propagatorType == 1 || propagatorType == 2 || propagatorType == 3 ) )
                {
                    // Progress
                    std::cout << std::endl << "Propagator: " << propagatorType + 1 << std::endl;

                    // Loop over integrators
                    integratorLimit = ( propagatorType == 7 ) ? 1 : 2; // only use the variable step size integrator
                    for ( int integratorType = 0; integratorType < integratorLimit; integratorType++ )
                    {
                        // Progress
                        std::cout << "Integrator: " << integratorType + 1 << std::endl;

                        // Loop over relativeTolerances and constant time step
                        valueLimit = ( integratorType == 0 ) ? relativeTolerances.size( ) : constantTimeStep.size( );
                        for ( int value = 0; value < valueLimit; value++ )
                        {
                            // Progress
                            std::cout << "Value: " << value + 1 << std::endl;

                            ///////////////////////     CREATE PROPAGATION SETTINGS         ////////////////////////////////////////////

                            // Propagator and integrator settings
                            boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings;
                            boost::shared_ptr< IntegratorSettings< > > integratorSettings;
                            boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                                    boost::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch );
                            if ( propagatorType == 7 )
                            {
                                // Propagator
                                translationalPropagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                                            centralBodies, accelerationModelMap, bodiesToPropagate, SatelliteInitialState, terminationSettings,
                                            cowell );

                                // Integrator
                                integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                            rungeKuttaVariableStepSize, simulationStartEpoch, 25.0,
                                            RungeKuttaCoefficients::rungeKuttaFehlberg78, 1e-5, 1e5,
                                            std::pow( 10, referenceTolerances ), std::pow( 10, referenceTolerances ) );//,
                                //                                    1, false, 0.8, 4.0, 0.1, &computeNewStepSize );
                            }
                            else
                            {
                                // Propagator
                                translationalPropagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                                            centralBodies, accelerationModelMap, bodiesToPropagate, SatelliteInitialState, simulationEndEpoch,
                                            static_cast< TranslationalPropagatorType >( propagatorType ) );

                                // Integrator
                                if ( integratorType == 0 )
                                {
                                    relativeTolerance = std::pow( 10, relativeTolerances.at( value ) );
                                    integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                                rungeKuttaVariableStepSize, simulationStartEpoch, 100.0,
                                                RungeKuttaCoefficients::rungeKuttaFehlberg56, 1e-5, 1e5,
                                                relativeTolerance, relativeTolerance );
//                                                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Constant( stateSizes.at( propagatorType ),
//                                                                                                      relativeTolerance ),
//                                                absoluteTolerances.at( propagatorType ) );//, 1, false, 0.8, 4.0, 0.1, &computeNewStepSize );
                                }
                                else if ( integratorType == 1 )
                                {
                                    integratorSettings = boost::make_shared< IntegratorSettings< > >(
                                                rungeKutta4, simulationStartEpoch, constantTimeStep.at( value ) );
                                }
                            }

                            // Add mass to propagator settings if 5th case
                            std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
                            propagatorSettingsVector.push_back( translationalPropagatorSettings );
                            if ( testCase == 5 )
                            {
                                // Create mass rate models
                                boost::shared_ptr< MassRateModelSettings > massRateModelSettings =
                                        boost::make_shared< FromThrustMassModelSettings >( true );
                                std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
                                massRateModels[ "Satellite" ] = createMassRateModel(
                                            "Satellite", massRateModelSettings, bodyMap, accelerationModelMap );

                                // Create settings for propagating the mass of the vehicle.
                                std::vector< std::string > bodiesWithMassToPropagate;
                                bodiesWithMassToPropagate.push_back( "Satellite" );

                                Eigen::VectorXd initialBodyMasses = Eigen::VectorXd( 1 );
                                initialBodyMasses( 0 ) = satelliteMass;

                                boost::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
                                        boost::make_shared< MassPropagatorSettings< double > >(
                                            bodiesWithMassToPropagate, massRateModels, initialBodyMasses, terminationSettings );

                                propagatorSettingsVector.push_back( massPropagatorSettings );
                            }
                            boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                                    boost::make_shared< MultiTypePropagatorSettings< double > >(
                                        propagatorSettingsVector, terminationSettings );

                            ///////////////////////     PROPAGATE ORBIT                     ////////////////////////////////////////////


                            // Simulate for 100 times to get a more accurate computation time
                            std::vector< double > computationTimes;
                            unsigned int numberOfSimulations = ( propagatorType == 7 ) ? 1 : 100;
                            boost::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator;
                            for ( unsigned int currentSimulation = 0; currentSimulation < numberOfSimulations; currentSimulation++ )
                            {
                                // Create simulation object and propagate dynamics
                                dynamicsSimulator = boost::make_shared< SingleArcDynamicsSimulator< > >(
                                            bodyMap, integratorSettings, propagatorSettings, true, false, false, false );
                                std::map< double, double > computationTimeMap = dynamicsSimulator->getCumulativeComputationTimeHistory( );
                                computationTimes.push_back( computationTimeMap.rbegin( )->second );
                            }
                            std::cout << "Average Computation Time: " << statistics::computeSampleMean( computationTimes ) << " s" << std::endl;

                            // Retrieve results
                            std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
                            std::map< double, Eigen::VectorXd > usmIntegrationResult;
                            if ( propagatorType > 3 && propagatorType < 7 )
                            {
                                usmIntegrationResult = dynamicsSimulator->getEquationsOfMotionNumericalSolutionRaw( );
                            }
                            std::map< double, unsigned int > functionEvaluationsMap = dynamicsSimulator->getCumulativeNumberOfFunctionEvaluations( );

                            ///////////////////////     PROVIDE OUTPUT TO FILES             ////////////////////////////////////////////

                            // Compute map of Kepler elements
                            Eigen::Vector6d cartesianState;
                            std::map< double, Eigen::VectorXd > keplerianIntegrationResult;
                            for ( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
                                  stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
                            {
                                // Retrieve current Cartesian state (convert to Earth-centered frame if needed)
                                cartesianState = stateIterator->second;
                                keplerianIntegrationResult[ stateIterator->first ] =
                                        convertCartesianToKeplerianElements( cartesianState, mainGravitationalParameter );
                            }

                            // Write perturbed satellite propagation history to file.
                            writeDataMapToTextFile( functionEvaluationsMap,
                                                    "eval" + nameAdditionPropagator[ propagatorType ] +
                                                    nameAdditionIntegrator[ integratorType ] +
                                                    "_" + std::to_string( simulationCase + 1 ) +
                                                    "_" + std::to_string( value + 1 ) + ".dat",
                                                    getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ) );

                            writeDataMapToTextFile( cartesianIntegrationResult,
                                                    "trajectory" + nameAdditionPropagator[ propagatorType ] +
                                                    nameAdditionIntegrator[ integratorType ] +
                                                    "_" + std::to_string( simulationCase + 1 ) +
                                                    "_" + std::to_string( value + 1 ) + ".dat",
                                                    getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ) );

                            writeDataMapToTextFile( keplerianIntegrationResult,
                                                    "orbit" + nameAdditionPropagator[ propagatorType ] +
                                                    nameAdditionIntegrator[ integratorType ] +
                                                    "_" + std::to_string( simulationCase + 1 ) +
                                                    "_" + std::to_string( value + 1 ) + ".dat",
                                                    getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ) );

                            if ( propagatorType > 3 && propagatorType < 7 )
                            {
                                writeDataMapToTextFile( usmIntegrationResult,
                                                        "usm" + nameAdditionPropagator[ propagatorType ] +
                                                        nameAdditionIntegrator[ integratorType ] +
                                                        "_" + std::to_string( simulationCase + 1 ) +
                                                        "_" + std::to_string( value + 1 ) + ".dat",
                                                        getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ) );
                            }

                            // Break loop if reference propagator
                            if ( propagatorType == 7 )
                                break;
                        }
                    }
                }
            }
        }

        // Break loop if only one case is to be run
        if ( singleTestCase )
            break;
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
