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
    double constantTimeStep;
    Eigen::Vector6d SatelliteInitialStateInKeplerianElements;

    // Select case
    //      0: Single aerobraking sweep
    //      1: Full aerobraking
    //      2: Circular orbit at LMO (Low Mars Orbit)
    //      3: Interplanetary trajectory
    int testCase = 1;
    std::vector< string > pathAdditionTestCase = { "aero", "aero_full", "circ", "inter" };
    switch ( testCase )
    {
    case 0: // Single aerobraking sweep
    {
        // Set simulation time settings.
        simulationDuration = 1.0 * physical_constants::JULIAN_DAY;
        constantTimeStep = 15.0;

        // Initial conditions
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 27185000;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.871988;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 158.7 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 23.4 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );
        break;
    }
    case 1: // Full aerobraking
    {
        // Set simulation time settings.
        simulationDuration = 150.0 * physical_constants::JULIAN_DAY;
        constantTimeStep = 150.0;

        // Initial conditions
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 27198500;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.871280;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 158.7 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 23.4 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );
        break;
    }
    case 2: // Circular orbit at LMO
    {
        // Set simulation time settings.
        simulationDuration = 27.5 * physical_constants::JULIAN_DAY;
        constantTimeStep = 25.0;

        // Initial conditions
        SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 3621000;
        SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.006904;
        SatelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
        SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 158.7 );
        SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 23.4 );
        SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 180.0 );
        break;
    }
    case 3: // Interplanetary trajectory
    {

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
    case 2:
    {
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( "Mars" );
        bodiesToCreate.push_back( "Earth" );
        break;
    }
    case 3:
    {
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( "Mercury" );
        bodiesToCreate.push_back( "Mars" );
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Jupiter" );
        bodiesToCreate.push_back( "Saturn" );
        bodiesToCreate.push_back( "Neptune" );
        bodiesToCreate.push_back( "Uranus" );
        break;
    }
    }

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

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    bodySettings[ "Mars" ]->gravityFieldSettings = boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    bodySettings[ "Mars" ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >( tabulatedAtmosphereFiles,
                                                                                                    atmosphereDependentVariables,
                                                                                                    atmosphereIndependentVariables,
                                                                                                    197.0, 1.3 );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Satellite" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Satellite" ]->setConstantBodyMass( 300.0 );

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
    double referenceAreaRadiation = 20.0;
    double radiationPressureCoefficient = 1.25;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
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
    bool keplerOrbit = false;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    if ( keplerOrbit )
    {
        accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    }
    else
    {
        accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 21, 21 ) );
        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            if ( i != 1 )
            {
                accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
            }
        }
        accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );

        accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    }

    // Add acceleration information
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Mars" );

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE INITIAL CONDITIONS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Convert to Cartesian elements
    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d SatelliteInitialState = convertKeplerianToCartesianElements(
                SatelliteInitialStateInKeplerianElements, marsGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             LOOP OVER PROPAGATORS                  ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Loop over propagators
    int integratorLimit;
    std::vector< string > nameAdditionPropagator = { "_cowell", "_encke", "_gauss", "_usm7", "_usm6", "_usmem", "_ref" };
    std::vector< string > nameAdditionIntegrator = { "_var", "_const" };
    for ( int propagatorType = 0; propagatorType < 7; propagatorType++ )
    {
        // Progress
        std::cout << "Propagator: " << propagatorType + 1 << std::endl;

        // Loop over integrators
        integratorLimit = ( propagatorType == 6 ) ? 1 : 2; // only use the variable step size integrator
        for ( int integratorType = 0; integratorType < integratorLimit; integratorType++ )
        {
            // Progress
            std::cout << "Integrator: " << integratorType + 1 << std::endl;

            ///////////////////////     CREATE PROPAGATION SETTINGS         ////////////////////////////////////////////

            // Propagator settings
            boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, SatelliteInitialState, simulationEndEpoch,
                      static_cast< TranslationalPropagatorType >( propagatorType ) );

            // Integrator settings
            boost::shared_ptr< IntegratorSettings< > > integratorSettings;
            if ( propagatorType == 6 )
            {
                integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                            rungeKuttaVariableStepSize, simulationStartEpoch, 25.0,
                            RungeKuttaCoefficients::rungeKuttaFehlberg78, 1e-3, 1e5, 1e-15, 1e-15 );
            }
            else
            {
                if ( integratorType == 0 )
                {
                    integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                rungeKuttaVariableStepSize, simulationStartEpoch, 100.0,
                                RungeKuttaCoefficients::rungeKuttaFehlberg56, 1e-1, 1e5, 1e-7, 1e-7 );
                }
                else if ( integratorType == 1 )
                {
                    integratorSettings = boost::make_shared< IntegratorSettings< > > (
                                rungeKutta4, simulationStartEpoch, constantTimeStep );
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
            std::map< double, double > functionEvaluationsMap = dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( );

            // Timing
            endTime = time( 0 );
            std::cout << "Propagation time: " << difftime( endTime, startTime ) << " s" << std::endl << std::endl;

            ///////////////////////     PROVIDE OUTPUT TO FILES             ////////////////////////////////////////////

            // Compute map of Kepler elements
            Eigen::Vector6d currentCartesianState;
            std::map< double, Eigen::VectorXd > keplerianIntegrationResult;
            std::map< double, Eigen::VectorXd > usmIntegrationResult;
            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
                 stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
            {
                // Retrieve current Cartesian state (convert to Earth-centered frame if needed)
                currentCartesianState = stateIterator->second;

                if ( centralBodies.at( 0 ) == "SSB" )
                {
                    currentCartesianState -=
                            bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( stateIterator->first );
                }

                keplerianIntegrationResult[ stateIterator->first ] =
                        convertCartesianToKeplerianElements( currentCartesianState, marsGravitationalParameter );

                if ( propagatorType == 3 )
                {
                    usmIntegrationResult[ stateIterator->first ] =
                            convertCartesianToUnifiedStateModelQuaternionsElements( currentCartesianState,
                                                                                    marsGravitationalParameter );
                }
                else if ( propagatorType == 4 )
                {
                    usmIntegrationResult[ stateIterator->first ] =
                            convertCartesianToUnifiedStateModelModifiedRodriguesParametersElements(
                                currentCartesianState, marsGravitationalParameter );
                }
                else if ( propagatorType == 5 )
                {
                    usmIntegrationResult[ stateIterator->first ] =
                            convertCartesianToUnifiedStateModelExponentialMapElements( currentCartesianState,
                                                                                       marsGravitationalParameter );
                }
            }

            // Write perturbed satellite propagation history to file.
            writeDataMapToTextFile( functionEvaluationsMap,
                                    "eval" + nameAdditionPropagator[ propagatorType ] +
                                    nameAdditionIntegrator[ integratorType ] + ".dat",
                                    getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ),
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );

            writeDataMapToTextFile( cartesianIntegrationResult,
                                    "trajectory" + nameAdditionPropagator[ propagatorType ] +
                                    nameAdditionIntegrator[ integratorType ] + ".dat",
                                    getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ),
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );

            writeDataMapToTextFile( keplerianIntegrationResult,
                                    "orbit" + nameAdditionPropagator[ propagatorType ] +
                                    nameAdditionIntegrator[ integratorType ] + ".dat",
                                    getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ),
                                    "",
                                    std::numeric_limits< double >::digits10,
                                    std::numeric_limits< double >::digits10,
                                    "," );

            if ( propagatorType > 2 && propagatorType < 6 )
            {
                writeDataMapToTextFile( usmIntegrationResult,
                                        "usm" + nameAdditionPropagator[ propagatorType ] +
                                        nameAdditionIntegrator[ integratorType ] + ".dat",
                                        getOutputPath( "Propagators/" + pathAdditionTestCase[ testCase ] ),
                                        "",
                                        std::numeric_limits< double >::digits10,
                                        std::numeric_limits< double >::digits10,
                                        "," );
            }
        }
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
