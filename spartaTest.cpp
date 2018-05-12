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

#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/Astrodynamics/Aerodynamics/rarefiedFlowAnalysis.h"

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/InputOutput/spartaInputOutput.h"

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "spartaTest.cpp" ).length( ) );
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
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::basic_mathematics;
    using namespace tudat::mathematical_constants;
    using namespace tudat::input_output;
    using namespace tudat::aerodynamics;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            TABULATED ATMOSPHERE          ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create a tabulated atmosphere object.
    std::string tabulatedAtmosphereFile = getAtmosphereTablesPath( ) + "MCDMeanAtmosphere.dat";
    std::vector< AtmosphereDependentVariables > dependentVariables = {
        density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere,
        gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere, molar_mass_dependent_atmosphere };
    TabulatedAtmosphere tabulatedAtmosphere = TabulatedAtmosphere( tabulatedAtmosphereFile,
                                                                   dependentVariables );

    // Extract atmosphere values
    std::map< double, Eigen::Vector3d > atmosphere;
    Eigen::Vector3d currentAtmosphere;
    int i = 0;
    for ( double h = 75; h <= 125; h++ )
    {
        currentAtmosphere[ 0 ] = tabulatedAtmosphere.getDensity( h * 1e3 );
        currentAtmosphere[ 1 ] = tabulatedAtmosphere.getPressure( h * 1e3 );
        currentAtmosphere[ 2 ] = tabulatedAtmosphere.getTemperature( h * 1e3 );
        atmosphere[ h * 1e3 ] = currentAtmosphere;
        i++;
    }

    writeDataMapToTextFile( atmosphere,
                            "atmosphere.dat",
                            getOutputPath( ),
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            INPUTS                        ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define inputs
    const std::string SPARTAExecutable = "'/Users/Michele/AE Software/SPARTA/src/spa_mac_mpi'";
    const int numberOfCores = 2;
    const std::string geometryFileUser = getSpartaDataPath( ) + "data/data.mro"; // check that it is not called data.shape
    const double referenceArea = 37.5;//3.12715;
    const double referenceLength = 1.0;
    const int referenceAxis = + 0; // axis opposite to free stream velocity and where reference aerodynamic area is taken
    Eigen::Vector3d momentReferencePoint = Eigen::Vector3d::Zero( );
    momentReferencePoint[ 1 ] = 1.25;
    const std::string simulationGases = "CO2"; // check that it matches gases in SPARTA atmosphere
    const double wallTemperature = 300;
    const double accomodationCoefficient = 1.0;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE INTERFACE              ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create analysis object.
    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultRarefiedFlowAltitudePoints( "Mars" );
    independentVariableDataPoints[ 1 ] = getDefaultRarefiedFlowMachPoints( "Full" );
    independentVariableDataPoints[ 2 ] = getDefaultRarefiedFlowAngleOfAttackPoints( );
//    independentVariableDataPoints[ 0 ] = { 100e3 };
//    independentVariableDataPoints[ 1 ] = { 17.1 };
//    independentVariableDataPoints[ 2 ] = { 0.0 * PI / 180.0 };

    // Generate database of aerodynamic coefficients.
    RarefiedFlowAnalysis coefficientInterface = RarefiedFlowAnalysis(
                SPARTAExecutable,
                independentVariableDataPoints,
                boost::make_shared< TabulatedAtmosphere >( tabulatedAtmosphere ),
                simulationGases,
                geometryFileUser,
                referenceArea,
                referenceLength,
                referenceAxis,
                momentReferencePoint,
                0.25,
                17.5,
                wallTemperature,
                accomodationCoefficient,
                false,
                "/usr/local/bin/mpirun",
                numberOfCores );

    // Write to files
    std::map< int, std::string > fileNamesMap;
    fileNamesMap[ 0 ] = input_output::getSpartaDataPath( ) + "coefficients/Cd.dat";
    fileNamesMap[ 1 ] = input_output::getSpartaDataPath( ) + "coefficients/Cs.dat";
    fileNamesMap[ 2 ] = input_output::getSpartaDataPath( ) + "coefficients/Cl.dat";
    fileNamesMap[ 3 ] = input_output::getSpartaDataPath( ) + "coefficients/Cm1.dat";
    fileNamesMap[ 4 ] = input_output::getSpartaDataPath( ) + "coefficients/Cm.dat";
    fileNamesMap[ 5 ] = input_output::getSpartaDataPath( ) + "coefficients/Cm3.dat";
    coefficientInterface.saveAerodynamicCoefficientsTables( fileNamesMap );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
