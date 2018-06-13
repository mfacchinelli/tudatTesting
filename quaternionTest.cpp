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
    const double simulationStartEpoch = 7.0 * tudat::physical_constants::JULIAN_YEAR +
            30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = 10.0 * tudat::physical_constants::JULIAN_DAY + simulationStartEpoch;

    // Define body settings for simulation
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );

    // Create body objects
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 1.0e3, simulationEndEpoch + 1.0e3 );

    // Reset ephemerides
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

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
    bodySettings[ "Mars" ]->gravityFieldSettings = boost::make_shared<
            FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    bodySettings[ "Mars" ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables,
                boundaryConditions, extrapolationValues );

    // Save body settings to map
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object
    bodyMap[ "Satellite" ] = boost::make_shared< Body >( );
    const double vehicleMass = 1000.0;
    bodyMap[ "Satellite" ]->setConstantBodyMass( vehicleMass );

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
    occultingBodies.push_back( "Mars" );
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
    accelerationsOfSatellite[ "Mars" ].push_back(
                boost::make_shared< SphericalHarmonicAccelerationSettings >( 21, 21 ) );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        if ( bodiesToCreate.at( i ) != "Mars" )
        {
            accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >(
                                                                              central_gravity ) );
        }
    }
    accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    // Add acceleration information
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Mars" );

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             INITIAL CONDITIONS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 23617.8637e3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.833437;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 62.8 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 280.0 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );

    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, marsGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             LOOP OVER NORMALIZATION METHODS        ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Loop over normalization methods
    //      0: no internal normalization
    //      1: normalization according to BOOK034 with rotational velocity
    //      2: normalization according to BOOK034 with constant factor
    //      3: normalization of derivative accoding to ART073
    //      4: normalization of quaternion and its derivative accoding to ART073
    unsigned int totalNumberOfMethods = 5;
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

        // Create simulation object and propagate dynamics
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false, true );
        std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > unifiedStateModelIntegrationResult =
                dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );
        std::map< double, unsigned int > functionEvaluationsMap = dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( );

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
