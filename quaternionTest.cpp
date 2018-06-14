/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Mathematics/Statistics/basicStatistics.h"

#include "Tudat/Astrodynamics/Propagators/quaternionNormalizationMethod.h"

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "quaternionTest.cpp" ).length( ) );
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

unsigned int tudat::propagators::QUATERNION_NORMALIZATION_METHOD;

//! Execute propagation of orbit of Asterix around Mars.
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
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings
    const double simulationStartEpoch = 7.0 * physical_constants::JULIAN_YEAR +
            30.0 * 6.0 * physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = 100.0 * tudat::physical_constants::JULIAN_DAY + simulationStartEpoch;

    // Define body settings for simulation
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 1.0e3, simulationEndEpoch + 1.0e3 );

    // Set ephemerides
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    // Generate gravitational and atmospheric body settings
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared<
            FromFileSphericalHarmonicsGravityFieldSettings >( ggm02s );
    bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< ExponentialAtmosphereSettings >(
                7050.0, 240.0, 1.225, 2.87e2 );

    // Save body settings to map
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object
    bodyMap[ "Satellite" ] = boost::make_shared< Body >( );
    const double satelliteMass = 1000.0;
    bodyMap[ "Satellite" ]->setConstantBodyMass( satelliteMass );

    // Aerodynamic coefficients from file
    std::map< int, std::string > aerodynamicCoefficientFiles;
    aerodynamicCoefficientFiles[ 0 ] = "/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/"
                                       "University/Master Thesis/Code/MATLAB/data/MRODragCoefficients.txt";
    aerodynamicCoefficientFiles[ 2 ] = "/Users/Michele/Library/Mobile Documents/com~apple~CloudDocs/"
                                       "University/Master Thesis/Code/MATLAB/data/MROLiftCoefficients.txt";

    // Create aerodynamic coefficient interface settings.
    const double referenceAreaAerodynamic = 37.5;
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = readTabulatedAerodynamicCoefficientsFromFiles(
                aerodynamicCoefficientFiles, referenceAreaAerodynamic,
                boost::assign::list_of( aerodynamics::angle_of_attack_dependent )( aerodynamics::altitude_dependent ),
                true, true );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = referenceAreaAerodynamic;
    double radiationPressureCoefficient = 1.25;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
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

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    accelerationsOfSatellite[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        if ( bodiesToCreate.at( i ) != "Earth" )
        {
            accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
    }
    accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfSatellite[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    // Add acceleration information
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Earth" );

    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             INITIAL CONDITIONS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d satelliteInitialStateInKeplerianElements;
    satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 26559.0e3;
    satelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.70;
    satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 63.2 );
    satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 206.3 );
    satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 281.6 );
    satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

    double marsGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements( satelliteInitialStateInKeplerianElements,
                                                                                     marsGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             LOOP OVER NORMALIZATION METHODS        ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Loop over normalization methods
    //      0: no internal normalization
    //      1: normalization according to BOOK034 with rotational velocity
    //      2: normalization according to BOOK034 with constant factor
    //      3: normalization of only quaternion accoding to ART073
    //      4: normalization of only derivative accoding to ART073
    //      5: normalization of both quaternion and derivative accoding to ART073
    unsigned int totalNumberOfMethods = 6;
    for ( unsigned int quaternionNormalizationMethod = 0;
          quaternionNormalizationMethod < ( totalNumberOfMethods + 1 ); quaternionNormalizationMethod++ )
    {
        // Set normalization method
        QUATERNION_NORMALIZATION_METHOD = quaternionNormalizationMethod;

        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////

        // Propagator and integrator settings
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
        boost::shared_ptr< IntegratorSettings< > > integratorSettings;
        if ( quaternionNormalizationMethod == totalNumberOfMethods )
        {
            // Propagator settings
            propagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch,
                      cowell );

            // Integrator settings
            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0, RungeKuttaCoefficients::rungeKuttaFehlberg78,
                        1e-5, 1e5, 1e-15, 1e-15 );
        }
        else
        {
            // Propagator settings
            propagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch,
                      unified_state_model_quaternions );

            // Integrator settings
            integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                        rungeKuttaVariableStepSize, simulationStartEpoch, 10.0, RungeKuttaCoefficients::rungeKuttaFehlberg56,
                        1e-5, 1e5, 1e-10, 1e-10 );
        }

        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////

        // Simulate for 10 times to get a more accurate computation time
        std::vector< double > computationTimes;
        unsigned int numberOfSimulations = ( quaternionNormalizationMethod != totalNumberOfMethods ) ? 10 : 1;
        boost::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator;
        for ( unsigned int currentSimulation = 0; currentSimulation < numberOfSimulations; currentSimulation++ )
        {
            // Create simulation object and propagate dynamics
            dynamicsSimulator = boost::make_shared< SingleArcDynamicsSimulator< > >(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false, true );
            std::map< double, double > computationTimeMap = dynamicsSimulator->getCumulativeComputationTimeHistory( );
            computationTimes.push_back( computationTimeMap.rbegin( )->second );
        }
        std::cout << "Average Computation Time: " << statistics::computeSampleMean( computationTimes ) << " s" << std::endl;

        // Retrieve results
        std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > unifiedStateModelIntegrationResult =
                dynamicsSimulator->getEquationsOfMotionNumericalSolutionRaw( );
        std::map< double, unsigned int > functionEvaluationsMap = dynamicsSimulator->getCumulativeNumberOfFunctionEvaluations( );

        ///////////////////////        PROVIDE OUTPUT TO FILES                       //////////////////////////////////////////

        // Write perturbed satellite propagation history in Cartesian elements to file
        input_output::writeDataMapToTextFile( cartesianIntegrationResult,
                                              "cartesianResult_" +
                                              std::to_string( quaternionNormalizationMethod ) + ".dat",
                                              getOutputPath( "Quaternions" ) );

        if ( quaternionNormalizationMethod != totalNumberOfMethods )
        {
            // Write perturbed satellite propagation history in USM7 elements to file
            input_output::writeDataMapToTextFile( unifiedStateModelIntegrationResult,
                                                  "usm7Result_" +
                                                  std::to_string( quaternionNormalizationMethod ) + ".dat",
                                                  getOutputPath( "Quaternions" ) );

            // Write function evaluation history to file
            input_output::writeDataMapToTextFile( functionEvaluationsMap,
                                                  "functionEvaluations_" +
                                                  std::to_string( quaternionNormalizationMethod ) + ".dat",
                                                  getOutputPath( "Quaternions" ) );
        }
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
