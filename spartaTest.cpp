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
        gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere };
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
//    independentVariableDataPoints[ 0 ] = getDefaultRarefiedFlowAltitudePoints( "Mars" );
//    independentVariableDataPoints[ 1 ] = getDefaultRarefiedFlowMachPoints( "High" );
//    independentVariableDataPoints[ 2 ] = getDefaultRarefiedFlowAngleOfAttackPoints( );
    independentVariableDataPoints[ 0 ] = { 100e3 };
    independentVariableDataPoints[ 1 ] = { 17.1 };
    independentVariableDataPoints[ 2 ] = { 0.0 * PI / 180.0 };

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
                accomodationCoefficient );

    // Write to files
    std::map< int, std::string > fileNamesMap;
    fileNamesMap[ 0 ] = input_output::getSpartaDataPath( ) + "coefficients/Cd.dat";
    fileNamesMap[ 1 ] = input_output::getSpartaDataPath( ) + "coefficients/Cs.dat";
    fileNamesMap[ 2 ] = input_output::getSpartaDataPath( ) + "coefficients/Cl.dat";
    fileNamesMap[ 3 ] = input_output::getSpartaDataPath( ) + "coefficients/Cm1.dat";
    fileNamesMap[ 4 ] = input_output::getSpartaDataPath( ) + "coefficients/Cm.dat";
    fileNamesMap[ 5 ] = input_output::getSpartaDataPath( ) + "coefficients/Cm3.dat";
    coefficientInterface.saveAerodynamicCoefficientsTables( fileNamesMap );

//    // Test basic properties of coefficient generator
//    std::cout << "Independent variables check: " << std::endl;
//    std::cout << coefficientInterface.getIndependentVariableNames( ).size( ) - 3 << std::endl;
//    std::cout << coefficientInterface.getIndependentVariableName( 0 ) << " " << altitude_dependent << std::endl;
//    std::cout << coefficientInterface.getIndependentVariableName( 1 ) << " " << mach_number_dependent << std::endl;
//    std::cout << coefficientInterface.getIndependentVariableName( 2 ) << " " << angle_of_attack_dependent << std::endl;

//    bool isVariableIndexTooHigh = 0;
//    try
//    {
//        coefficientInterface.getIndependentVariableName( 3 );
//    }
//    catch ( std::runtime_error )
//    {
//        isVariableIndexTooHigh = 1;
//    }
//    std::cout << isVariableIndexTooHigh << std::endl << std::endl;

//    // Allocate memory for independent variables to pass to analysis for retrieval.
//    boost::array< int, 3 > independentVariables;
//    independentVariables[ 0 ] = 0;
//    independentVariables[ 1 ] = 0;
//    independentVariables[ 2 ] = 0;
//    std::vector< double > independentVariablesVector( 3 );
//    std::vector< double > interpolatingIndependentVariablesVector( 3 );
//    const double expectedValueOfForceCoefficient = 0.0;

//    // Declare local test variables.
//    Eigen::Vector6d aerodynamicCoefficients_ = Eigen::Vector6d::Zero( );
//    double forceCoefficient_;

//    // Iterate over all angles of attack to verify sphere coefficients. Total force coefficient
//    // should be one; all moment coefficients should be zero.
//    // The functionality is tested directly from the generator, as well as from the
//    // coefficient interface, both interpolated at the nodes, and halfway between the nodes.
//    for ( int i = 0; i < coefficientInterface.getNumberOfValuesOfIndependentVariable( 0 ); i++ )
//    {
//        independentVariables[ 0 ] = i;
//        independentVariablesVector[ 0 ] = coefficientInterface.getIndependentVariablePoint( 0, i );
//        if ( i < coefficientInterface.getNumberOfValuesOfIndependentVariable( 0 ) - 1 )
//        {
//            interpolatingIndependentVariablesVector[ 0 ] =
//                    coefficientInterface.getIndependentVariablePoint( 0, i ) + 0.5 * (
//                        coefficientInterface.getIndependentVariablePoint( 0, i + 1 ) -
//                        coefficientInterface.getIndependentVariablePoint( 0, i ) );
//        }

//        for ( int j = 0; j <
//              coefficientInterface.getNumberOfValuesOfIndependentVariable( 1 ); j++ )
//        {
//            independentVariables[ 1 ] = j;
//            independentVariablesVector[ 1 ] =
//                    coefficientInterface.getIndependentVariablePoint( 1, j );
//            if ( j < coefficientInterface.getNumberOfValuesOfIndependentVariable( 1 ) - 1 )
//            {
//                interpolatingIndependentVariablesVector[ 1 ] =
//                        coefficientInterface.getIndependentVariablePoint( 1, j ) + 0.5 * (
//                            coefficientInterface.getIndependentVariablePoint( 1, j + 1 ) -
//                            coefficientInterface.getIndependentVariablePoint( 1, j ) );
//            }

//            for ( int k = 0; k <
//                  coefficientInterface.getNumberOfValuesOfIndependentVariable( 2 ); k++ )
//            {
//                independentVariables[ 2 ] = k;
//                independentVariablesVector[ 2 ] =
//                        coefficientInterface.getIndependentVariablePoint( 2, k );
//                if ( k < coefficientInterface.getNumberOfValuesOfIndependentVariable( 2 ) - 1 )
//                {
//                    interpolatingIndependentVariablesVector[ 2 ] =
//                            coefficientInterface.getIndependentVariablePoint( 2, k ) + 0.5 * (
//                                coefficientInterface.getIndependentVariablePoint( 2, k + 1 ) -
//                                coefficientInterface.getIndependentVariablePoint( 2, k ) );
//                }
//                std::cout << std::endl << "I: " << i << " J: " << j << " K: " << k << std::endl;

//                // Retrieve aerodynamic coefficients.
//                aerodynamicCoefficients_ =
//                        coefficientInterface.getAerodynamicCoefficientsDataPoint(
//                            independentVariables );
//                forceCoefficient_ = ( aerodynamicCoefficients_.head( 3 ) ).norm( );

//                // Test if the computed force coefficient corresponds to the expected value
//                // within the specified tolerance.
//                std::cout << std::endl << "1: " << forceCoefficient_ - expectedValueOfForceCoefficient << std::endl;

//                // Test if the computed moment coefficients correspond to the expected value (0.0)
//                // within the specified tolerance.
//                std::cout << "2: " << aerodynamicCoefficients_( 3 ) << std::endl;
//                std::cout << "3: " << aerodynamicCoefficients_( 4 ) << std::endl;
//                std::cout << "4: " << aerodynamicCoefficients_( 5 ) << std::endl;

//                // Retrieve aerodynamic coefficients from coefficient interface.
//                coefficientInterface.updateCurrentCoefficients( independentVariablesVector );

//                aerodynamicCoefficients_ =
//                        coefficientInterface.getCurrentAerodynamicCoefficients( );
//                forceCoefficient_ = ( aerodynamicCoefficients_.head( 3 ) ).norm( );

//                // Test if the computed force coefficient corresponds to the expected value
//                // within the specified tolerance.
//                std::cout << std::endl << "5: " << forceCoefficient_ - expectedValueOfForceCoefficient << std::endl;

//                // Test if the computed moment coefficients correspond to the expected value (0.0)
//                // within the specified tolerance.
//                std::cout << "6: " << aerodynamicCoefficients_( 3 ) << std::endl;
//                std::cout << "7: " << aerodynamicCoefficients_( 4 ) << std::endl;
//                std::cout << "8: " << aerodynamicCoefficients_( 5 ) << std::endl;

//                // Retrieve aerodynamic coefficients from coefficient interface.
//                coefficientInterface.updateCurrentCoefficients(
//                            interpolatingIndependentVariablesVector );

//                aerodynamicCoefficients_ =
//                        coefficientInterface.getCurrentAerodynamicCoefficients( );
//                forceCoefficient_ = ( aerodynamicCoefficients_.head( 3 ) ).norm( );

//                // Test if the computed force coefficient corresponds to the expected value
//                // within the specified tolerance.
//                std::cout << std::endl << "9: " << forceCoefficient_ - expectedValueOfForceCoefficient << std::endl;

//                // Test if the computed moment coefficients correspond to the expected value (0.0)
//                // within the specified tolerance.
//                std::cout << "10: " << aerodynamicCoefficients_( 3 ) << std::endl;
//                std::cout << "11: " << aerodynamicCoefficients_( 4 ) << std::endl;
//                std::cout << "12: " << aerodynamicCoefficients_( 5 ) << std::endl;
//            }
//        }
//    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
