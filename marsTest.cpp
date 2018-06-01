/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "marsTest.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = 7.0 * tudat::physical_constants::JULIAN_YEAR +
            30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = 1.0 * tudat::physical_constants::JULIAN_DAY + simulationStartEpoch;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Saturn" );

    // Atmosphere
    std::string atmosphereFile = "/Users/Michele/GitHub/tudat/tudatBundle/tudat/Tudat/External/AtmosphereTables/MCDMeanAtmosphere.dat";
    std::vector< aerodynamics::AtmosphereDependentVariables > atmosphereDependentVariables;
    for ( unsigned int i = 0; i < 6; i++ )
    {
        atmosphereDependentVariables.push_back( static_cast< aerodynamics::AtmosphereDependentVariables >( i ) );
    }

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    bodySettings[ "Mars" ]->gravityFieldSettings =
            boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    bodySettings[ "Mars" ]->atmosphereSettings =
            boost::make_shared< TabulatedAtmosphereSettings >( atmosphereFile, atmosphereDependentVariables );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Satellite" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Satellite" ]->setConstantBodyMass( 1000.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 37.5;
    double aerodynamicCoefficient = 2.2;
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 37.5;
    double radiationPressureCoefficient = 1.0;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Satellite" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Satellite", bodyMap ) );

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
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        if ( bodiesToCreate.at( i ) != "Mars" )
        {
            accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back(
                        boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
        }
    }
    accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );
    accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                      basic_astrodynamics::aerodynamic ) );

    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Mars" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 23477500;//25945000;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.848152;//0.865099;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 93.0 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 158.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 180.0 );

    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, marsGravitationalParameter );

    // Propagator settings
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch,
              unified_state_model_exponential_map ); // unified_state_model_exponential_map

    // Integrator settings
//    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
//            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                rungeKuttaVariableStepSize, simulationStartEpoch, 1.0,
//                RungeKuttaCoefficients::rungeKuttaFehlberg56, 1e-3, 1e3, 1e-12, 1e-12 );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > > (
                rungeKutta4, simulationStartEpoch, 1.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > rawIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );
    std::map< double, Eigen::VectorXd > keplerianIntegrationResult;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute map of Kepler elements
    Eigen::Vector6d currentCartesianState;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
         stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
    {
        // Retrieve current Cartesian state (convert to Mars-centered frame if needed)
        currentCartesianState = stateIterator->second;
        keplerianIntegrationResult[ stateIterator->first ] =
                convertCartesianToKeplerianElements( currentCartesianState, marsGravitationalParameter );
    }

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( cartesianIntegrationResult,
                                          "cart.dat", getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( rawIntegrationResult,
                                          "raw.dat", getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

