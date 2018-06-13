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

#include "SatellitePropagatorExamples/applicationOutput.h"

//! Execute propagation of orbit of Satellite around the Earth.
int main()
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
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 10*tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
//    bodiesToCreate.push_back( "Moon" );
//    bodiesToCreate.push_back( "Mars" );
//    bodiesToCreate.push_back( "Venus" );
//    bodiesToCreate.push_back( "Jupiter" );
//    bodiesToCreate.push_back( "Saturn" );
//    bodiesToCreate.push_back( "Neptune" );
//    bodiesToCreate.push_back( "Uranus" );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 ); // NRLMSISE-00
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Satellite" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Satellite" ]->setConstantBodyMass( 300.0 );

    // Aerodynamic parameters
    const double referenceAreaAerodynamic = 37.5;
    Eigen::Vector3d aerodynamicCoefficient = Eigen::Vector3d::Zero( );
    aerodynamicCoefficient( 0 ) = 2.2; // drag coefficient

    // Create aerodynamic coefficient interface settings.
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceAreaAerodynamic, aerodynamicCoefficient, 1, 1 );

    bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 20;
    double radiationPressureCoefficient = 1.25;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
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
    accelerationsOfSatellite[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 21, 21 ) );

    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >(
                                                                         basic_astrodynamics::central_gravity ) );
    }
    accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationsOfSatellite[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::aerodynamic ) );

    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Earth" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Satellite.
    Eigen::Vector6d SatelliteInitialStateInKeplerianElements;
    SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = -10763333.333333;
    SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 1.6;
    SatelliteInitialStateInKeplerianElements( inclinationIndex ) = 0.7853981634;
    SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = 0;
    SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = 0;
    SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = 5.387521;

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d SatelliteInitialState = convertKeplerianToCartesianElements(
                SatelliteInitialStateInKeplerianElements, earthGravitationalParameter );

    // Propagation termination conditions.
    boost::shared_ptr< PropagationTerminationSettings > terminationSettingsTime = boost::make_shared<
            PropagationTimeTerminationSettings >( simulationEndEpoch );
    boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettingAltitude =
            boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Satellite", "Earth" );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettingsAltitude = boost::make_shared<
            PropagationDependentVariableTerminationSettings >( dependentVariableSettingAltitude, 10.0E+03, true );
    std::vector< boost::shared_ptr< PropagationTerminationSettings > > propagationTerminationSettings;
    propagationTerminationSettings.push_back( terminationSettingsTime );
    propagationTerminationSettings.push_back( terminationSettingsAltitude );

    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationHybridTerminationSettings >( propagationTerminationSettings, true );

    // Dependent variables
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable, "Satellite", "Earth" ) );
    dependentVariables.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable, "Satellite", "Earth" ) );
    dependentVariables.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::aerodynamic, "Satellite", "Earth", false ) );

    // Propagator settings
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, SatelliteInitialState, terminationSettings,
              unified_state_model_quaternions, boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

    // Integrator settings
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                rungeKuttaVariableStepSize, simulationStartEpoch, 1.0,
                RungeKuttaCoefficients::rungeKuttaFehlberg56, 0.005, 10.0, 1e-15, 1e-15 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > keplerianIntegrationResult;
    std::map< double, Eigen::VectorXd > dependentVariableSolution = dynamicsSimulator.getDependentVariableHistory( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute map of Kepler elements
    Eigen::Vector6d currentCartesianState;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
         stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
    {
        // Retrieve current Cartesian state (convert to Earth-centered frame if needed)
        currentCartesianState = stateIterator->second;

        if( centralBodies.at( 0 ) == "SSB" )
        {
            currentCartesianState -=
                    bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( stateIterator->first );
        }

        keplerianIntegrationResult[ stateIterator->first ] =
                convertCartesianToKeplerianElements( currentCartesianState, earthGravitationalParameter );
    }

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( cartesianIntegrationResult,
                                          "trajectory.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( keplerianIntegrationResult,
                                          "orbit.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( dependentVariableSolution,
                                          "dependent.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

