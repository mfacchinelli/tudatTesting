/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <ctime>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

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
    ///////////////////////            DEFINE TEST CASES             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Predefine variables
    double simulationDuration = 0.0;
    std::vector< double > constantTimeStep;
    std::string simulationCentralBody;
    int limitingSphericalHarminics;
    Eigen::Vector6d SatelliteInitialStateInKeplerianElements;

    // Variables that can be optionally overwritten
    std::vector< int > tolerances = { -13, -12, -11, -10, -9, -8, -7 };
    int referenceTolerance = -15;

    // Select case
    //      0: Aerocapture
    //      1: Full aerobraking
    //      2: Interplanetary trajectory
    //      3: Circular orbit at LEO (Low Earth Orbit)
    //      4: Molniya orbit
    int testCase = 5;
    std::vector< std::string > pathAdditionTestCase = { "aero", "aero_full", "inter", "circ", "moln" };
    switch ( testCase )
    {
    case 0: // Aerocapture
    {
        // Set simulation time settings
        simulationDuration = 0.625 * physical_constants::JULIAN_DAY;
        constantTimeStep = { 0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0 };

        // Set simulation central body
        simulationCentralBody = "Mars";
        limitingSphericalHarminics = 21;

        // Initial conditions
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = -17305000.000000;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 1.2;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 23.4 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 230.0 );
        break;
    }
    case 1: // Full aerobraking
    {
        // Set simulation time settings.
        simulationDuration = 100.0 * physical_constants::JULIAN_DAY;
        constantTimeStep = { 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0 };

        // Set simulation central body
        simulationCentralBody = "Mars";
        limitingSphericalHarminics = 21;

        // Initial conditions
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 27221000;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.869733;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = 0.0;//convertDegreesToRadians( 93.0 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;//convertDegreesToRadians( 158.7 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;//convertDegreesToRadians( 23.4 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = 0.0;//convertDegreesToRadians( 180.0 );
        break;
    }
    case 2: // Interplanetary trajectory
    {
        // Set simulation time settings.
        simulationDuration = 50.0 * physical_constants::JULIAN_DAY;
        constantTimeStep = { 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0 };

        // Set simulation central body
        simulationCentralBody = "Sun";

        // Use less stringent tolerances
        tolerances = { -11, -10, -9, -8, -7, -6, -5 };
        referenceTolerance = -13;

        // Initial conditions
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 845508000;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.059136;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 2.5 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 328.6 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 147.8 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
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
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6936000.0;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.0;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 28.5 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 194.8 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 272.3 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
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
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 26559000.0;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.70;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 63.2 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 206.3 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 281.6 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
        break;
    }
    case 5: // Vittaldev low-thrust orbit
    {
        // Set simulation time settings.
        simulationDuration = 10.0 * physical_constants::JULIAN_DAY;
        constantTimeStep = { 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0 };

        // Set simulation central body
        simulationCentralBody = "Earth";
        limitingSphericalHarminics = 21;

        // Initial conditions
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6378.1363e3 + 1000e3;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.0;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 0.0 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 0.0 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 0.0 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
        break;
    }
    default:
    {
        std::cerr << "Mode not recognized." << std::endl;
    }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = 7.0 * physical_constants::JULIAN_YEAR +
            30.0 * 6.0 * physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = simulationDuration + simulationStartEpoch;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    switch ( testCase )
    {
    case 0:
    case 1:
    {
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( simulationCentralBody );
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
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( simulationCentralBody );
        bodiesToCreate.push_back( "Moon" );
        break;
    }
    }

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
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
        tabulatedAtmosphereFiles[ 0 ] = getAtmosphereTablesPath( ) +
                "MCDMeanAtmosphereTimeAverage/density.dat";
        tabulatedAtmosphereFiles[ 1 ] = getAtmosphereTablesPath( ) +
                "MCDMeanAtmosphereTimeAverage/pressure.dat";
        tabulatedAtmosphereFiles[ 2 ] = getAtmosphereTablesPath( ) +
                "MCDMeanAtmosphereTimeAverage/temperature.dat";
        tabulatedAtmosphereFiles[ 3 ] = getAtmosphereTablesPath( ) +
                "MCDMeanAtmosphereTimeAverage/gasConstant.dat";
        tabulatedAtmosphereFiles[ 4 ] = getAtmosphereTablesPath( ) +
                "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";
        std::vector< AtmosphereDependentVariables > atmosphereDependentVariables = {
            density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere,
            gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere };
        std::vector< AtmosphereIndependentVariables > atmosphereIndependentVariables = {
            longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };
        std::vector< interpolators::BoundaryInterpolationType > boundaryConditions = {
            interpolators::use_boundary_value, interpolators::use_boundary_value, interpolators::use_default_value };
        std::vector< double > extrapolationValues = { 0.0, 0.0, 186.813, 8183.0, 1.667 };

        bodySettings[ simulationCentralBody ]->gravityFieldSettings = boost::make_shared<
                FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
        bodySettings[ simulationCentralBody ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                    tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables,
                    boundaryConditions, extrapolationValues );
        break;
    }
    case 2:
    {
        break;
    }
    case 3:
    case 4:
    case 5:
    {
        bodySettings[ simulationCentralBody ]->gravityFieldSettings = boost::make_shared<
                FromFileSphericalHarmonicsGravityFieldSettings >( ggm02s );
        bodySettings[ simulationCentralBody ]->atmosphereSettings = boost::make_shared< ExponentialAtmosphereSettings >(
                    7050, 240.0, 1.225, 2.87e2 );
//        bodySettings[ simulationCentralBody ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
//                    getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );
        break;
    }
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Satellite" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Satellite" ]->setConstantBodyMass( 1000.0 );

    // Aerodynamic parameters
    const double referenceAreaAerodynamic = 37.5;

    // Aerodynamic coefficients from file
    std::map< int, std::string > aerodynamicCoefficientFiles;
    aerodynamicCoefficientFiles[ 0 ] = "/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/"
                                       "University/Master Thesis/Code/MATLAB/data/MRODragCoefficients.txt";
    aerodynamicCoefficientFiles[ 2 ] = "/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/"
                                       "University/Master Thesis/Code/MATLAB/data/MROLiftCoefficients.txt";

    // Create aerodynamic coefficient interface settings.
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
                aerodynamicCoefficientFiles, referenceAreaAerodynamic,
                boost::assign::list_of( aerodynamics::angle_of_attack_dependent )( aerodynamics::altitude_dependent ),
                true, true );

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
    bodyMap[ "Satellite" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    SatelliteRadiationPressureSettings, "Satellite", bodyMap ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    switch ( testCase )
    {
    case 0:
    case 1:
    case 3:
    case 4:
    case 5:
    {
        bool keplerOrbit = false;
        if ( keplerOrbit )
        {
            accelerationsOfSatellite[ simulationCentralBody ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
        else
        {
            accelerationsOfSatellite[ simulationCentralBody ].push_back(
                        boost::make_shared< SphericalHarmonicAccelerationSettings >( limitingSphericalHarminics,
                                                                                     limitingSphericalHarminics ) );
            for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
            {
                if ( bodiesToCreate.at( i ) != simulationCentralBody )
                {
                    accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
                }
            }
            accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );

            accelerationsOfSatellite[ simulationCentralBody ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
        }

        // Add thrust acceleration if 5th case
        if ( testCase == 5 )
        {
            // Define thrust settings
            double thrustMagnitude = 9.81;
            double specificImpulse = 5000.0;
            boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
                    boost::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                        simulationCentralBody, true, false );
            boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
                    boost::make_shared< ConstantThrustEngineSettings >(
                        thrustMagnitude, specificImpulse );

            // Define thrust acceleration settings
            accelerationsOfSatellite[ "Satellite" ].push_back(
                        boost::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings,
                                                                          thrustMagnitudeSettings) );
        }
        break;
    }
    case 2:
    {
        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
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
                SatelliteInitialStateInKeplerianElements, mainGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             LOOP OVER PROPAGATORS                  ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Loop over propagators
    int integratorLimit, valueLimit;
    double tolerance;
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

                // Loop over tolerances and constant time step
                valueLimit = ( integratorType == 0 ) ? tolerances.size( ) : constantTimeStep.size( );
                for ( int value = 0; value < valueLimit; value++ )
                {
                    // Progress
                    std::cout << "Value: " << value + 1 << std::endl;

                    ///////////////////////     CREATE PROPAGATION SETTINGS         ////////////////////////////////////////////

                    // Propagator and integrator settings
                    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
                    boost::shared_ptr< IntegratorSettings< > > integratorSettings;
                    if ( propagatorType == 7 )
                    {
                        // Propagator
                        propagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >
                                ( centralBodies, accelerationModelMap, bodiesToPropagate, SatelliteInitialState, simulationEndEpoch,
                                  cowell );

                        // Integrator
                        integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                    rungeKuttaVariableStepSize, simulationStartEpoch, 25.0,
                                    RungeKuttaCoefficients::rungeKuttaFehlberg78, 1e-5, 1e5,
                                    std::pow( 10, referenceTolerance ), std::pow( 10, referenceTolerance ) );
                    }
                    else
                    {
                        // Propagator
                        propagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >
                                ( centralBodies, accelerationModelMap, bodiesToPropagate, SatelliteInitialState, simulationEndEpoch,
                                  static_cast< TranslationalPropagatorType >( propagatorType ) );

                        // Integrator
                        if ( integratorType == 0 )
                        {
                            tolerance = std::pow( 10, tolerances.at( value ) );
                            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        rungeKuttaVariableStepSize, simulationStartEpoch, 100.0,
                                        RungeKuttaCoefficients::rungeKuttaFehlberg56, 1e-5, 1e5,
                                        tolerance, tolerance );
                        }
                        else if ( integratorType == 1 )
                        {
                            integratorSettings = boost::make_shared< IntegratorSettings< > > (
                                        rungeKutta4, simulationStartEpoch, constantTimeStep.at( value ) );
                        }
                    }

                    ///////////////////////     PROPAGATE ORBIT                     ////////////////////////////////////////////

                    // Timeing
                    time_t startTime, endTime;
                    startTime = time( 0 );

                    // Create simulation object and propagate dynamics.
                    SingleArcDynamicsSimulator< > dynamicsSimulator(
                                bodyMap, integratorSettings, propagatorSettings, true, false, false );
                    std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                    std::map< double, Eigen::VectorXd > usmIntegrationResult;
                    if ( propagatorType > 3 && propagatorType < 7 )
                    {
                        usmIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );
                    }
                    std::map< double, double > functionEvaluationsMap = dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( );

                    // Timing
                    endTime = time( 0 );
                    std::cout << "Propagation time: " << difftime( endTime, startTime ) << " s" << std::endl;

                    ///////////////////////     PROVIDE OUTPUT TO FILES             ////////////////////////////////////////////

                    // Compute map of Kepler elements
                    Eigen::Vector6d cartesianState;
                    std::map< double, Eigen::VectorXd > keplerianIntegrationResult;
                    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
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
                                            "_" + std::to_string( value + 1 ) + ".dat",
                                            getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ),
                                            "",
                                            std::numeric_limits< double >::digits10,
                                            std::numeric_limits< double >::digits10,
                                            "," );

                    writeDataMapToTextFile( cartesianIntegrationResult,
                                            "trajectory" + nameAdditionPropagator[ propagatorType ] +
                                            nameAdditionIntegrator[ integratorType ] +
                                            "_" + std::to_string( value + 1 ) + ".dat",
                                            getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ),
                                            "",
                                            std::numeric_limits< double >::digits10,
                                            std::numeric_limits< double >::digits10,
                                            "," );

                    writeDataMapToTextFile( keplerianIntegrationResult,
                                            "orbit" + nameAdditionPropagator[ propagatorType ] +
                                            nameAdditionIntegrator[ integratorType ] +
                                            "_" + std::to_string( value + 1 ) + ".dat",
                                            getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ),
                                            "",
                                            std::numeric_limits< double >::digits10,
                                            std::numeric_limits< double >::digits10,
                                            "," );

                    if ( propagatorType > 3 && propagatorType < 7 )
                    {
                        writeDataMapToTextFile( usmIntegrationResult,
                                                "usm" + nameAdditionPropagator[ propagatorType ] +
                                                nameAdditionIntegrator[ integratorType ] +
                                                "_" + std::to_string( value + 1 ) + ".dat",
                                                getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ),
                                                "",
                                                std::numeric_limits< double >::digits10,
                                                std::numeric_limits< double >::digits10,
                                                "," );
                    }

                    // Break loop if reference propagator
                    if ( propagatorType == 7 )
                    {
                        break;
                    }
                }
            }
        }
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
