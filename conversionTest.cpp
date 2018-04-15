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

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/mapTextFileReader.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "conversionTest.cpp" ).length( ) );
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

//! Execute propagation of orbit of Satellite around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CREATE KEPLER ELEMENT STATE        //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create list of Keplerian angles
    std::vector< double > inclinationAngles;
    std::vector< double > rightAscensionOfAscendingNodeAngles;
    std::vector< double > argumentOfPeriapsisAngles;
    std::vector< double > trueAnomalyAngles;
    double angle = 0.0;
    while ( angle < 360.0 )
    {
        // Add inclination to list if below definition boundary
        if ( angle < 180.0 )
        {
            inclinationAngles.push_back( unit_conversions::convertDegreesToRadians( angle ) );
        }
        else if ( angle == 180.0 )
        {
            inclinationAngles.push_back( unit_conversions::convertDegreesToRadians( angle ) - 0.001 );
        }

        // Add other angles up to 360 degrees
        rightAscensionOfAscendingNodeAngles.push_back( unit_conversions::convertDegreesToRadians( angle ) );
        argumentOfPeriapsisAngles.push_back( unit_conversions::convertDegreesToRadians( angle ) );
        trueAnomalyAngles.push_back( unit_conversions::convertDegreesToRadians( angle ) );

        // Next step
        angle += 7.5;
    }

    // Create Keplerian elements vector
    Eigen::Vector6d keplerianElements = Eigen::Vector6d::Zero( );
    keplerianElements( semiMajorAxisIndex ) = 1.5e11;
    keplerianElements( eccentricityIndex ) = 0.25;

    // Standard gravitational parameter
    const double centralBodyGravitationalParameter = 1.32712440018e20; // [m^3/s^2]

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CONVERT TO AND FROM USMEM          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Timeing
    time_t tstart, tend;
    tstart = time( 0 );

    // Create map of results
    std::map< int, Eigen::Vector6d > keplerianElementsInput;
    std::map< int, Eigen::Vector6d > keplerianElementsUSM7Output;
    std::map< int, Eigen::Vector6d > keplerianElementsUSM6Output;
    std::map< int, Eigen::Vector6d > keplerianElementsUSMEMOutput;
    std::map< int, Eigen::Vector7d > unifiedStateModelQuaternionsElementsOutput;
    std::map< int, Eigen::Vector7d > unifiedStateModelModifiedRodriguesParametersElementsOutput;
    std::map< int, Eigen::Vector6d > unifiedStateModelExponentialMapElementsOutput;

    // Loop over angles
    int loopIndex = 1;
    double tolerance = 1e-15;
    Eigen::Vector6d convertedKeplerianElementsUSM7 = Eigen::Vector6d::Zero( );
    Eigen::Vector6d convertedKeplerianElementsUSM6 = Eigen::Vector6d::Zero( );
    Eigen::Vector6d convertedKeplerianElementsUSMEM = Eigen::Vector6d::Zero( );
    Eigen::Vector7d convertedUnifiedStateModelQuaternionsElements = Eigen::Vector7d::Zero( );
    Eigen::Vector7d convertedUnifiedStateModelModifiedRodriguesParametersElements = Eigen::Vector7d::Zero( );
    Eigen::Vector6d convertedUnifiedStateModelExponentialMapElements = Eigen::Vector6d::Zero( );
    for ( std::vector< double >::const_iterator i = inclinationAngles.begin( ); i != inclinationAngles.end( ); ++i )
    {
        for ( std::vector< double >::const_iterator O = rightAscensionOfAscendingNodeAngles.begin( ); O != rightAscensionOfAscendingNodeAngles.end( ); ++O )
        {
            for ( std::vector< double >::const_iterator o = argumentOfPeriapsisAngles.begin( ); o != argumentOfPeriapsisAngles.end( ); ++o )
            {
                for ( std::vector< double >::const_iterator t = trueAnomalyAngles.begin( ); t != trueAnomalyAngles.end( ); ++t )
                {
                    // Complete Keplerian elements vector
                    keplerianElements( inclinationIndex ) = *i;
                    keplerianElements( longitudeOfAscendingNodeIndex ) = *O;
                    keplerianElements( argumentOfPeriapsisIndex ) = *o;
                    keplerianElements( trueAnomalyIndex ) = *t;

                    // Check for inconsistencies
                    if ( std::abs( keplerianElements( inclinationIndex ) ) < tolerance )
                    {
                        keplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
                    }

                    // Convert to USM7 and back
                    Eigen::Vector6d cartesianElements = convertKeplerianToCartesianElements(
                                keplerianElements, centralBodyGravitationalParameter );
                    convertedUnifiedStateModelQuaternionsElements = convertCartesianToUnifiedStateModelQuaternionsElements(
                                cartesianElements, centralBodyGravitationalParameter );

                    Eigen::Vector6d convertedCartesianElementsUSM7 = convertUnifiedStateModelQuaternionsToCartesianElements(
                                convertedUnifiedStateModelQuaternionsElements, centralBodyGravitationalParameter );
                    convertedKeplerianElementsUSM7 = convertCartesianToKeplerianElements(
                                convertedCartesianElementsUSM7, centralBodyGravitationalParameter );

//                    // Convert to USM7 and back
//                    convertedUnifiedStateModelQuaternionsElements = convertKeplerianToUnifiedStateModelQuaternionsElements(
//                                keplerianElements, centralBodyGravitationalParameter );

//                    convertedKeplerianElementsUSM7 = convertUnifiedStateModelQuaternionsToKeplerianElements(
//                                convertedUnifiedStateModelQuaternionsElements, centralBodyGravitationalParameter );

                    // Convert to USM6 and back
                    convertedUnifiedStateModelModifiedRodriguesParametersElements =
                            convertKeplerianToUnifiedStateModelModifiedRodriguesParametersElements(
                                keplerianElements, centralBodyGravitationalParameter );

                    convertedKeplerianElementsUSM6 =
                            convertUnifiedStateModelModifiedRodriguesParametersToKeplerianElements(
                                convertedUnifiedStateModelModifiedRodriguesParametersElements,
                                centralBodyGravitationalParameter );

                    // Convert to USMEM and back
                    convertedUnifiedStateModelExponentialMapElements = convertKeplerianToUnifiedStateModelExponentialMapElements(
                                keplerianElements, centralBodyGravitationalParameter );

                    convertedKeplerianElementsUSMEM = convertUnifiedStateModelExponentialMapToKeplerianElements(
                                convertedUnifiedStateModelExponentialMapElements, centralBodyGravitationalParameter );

                    // Save to map
                    keplerianElementsInput[ loopIndex ] = keplerianElements;
                    keplerianElementsUSM7Output[ loopIndex ] = convertedKeplerianElementsUSM7;
                    keplerianElementsUSM6Output[ loopIndex ] = convertedKeplerianElementsUSM6;
                    keplerianElementsUSMEMOutput[ loopIndex ] = convertedKeplerianElementsUSMEM;
                    unifiedStateModelQuaternionsElementsOutput[ loopIndex ] = convertedUnifiedStateModelQuaternionsElements;
                    unifiedStateModelModifiedRodriguesParametersElementsOutput[ loopIndex ] = convertedUnifiedStateModelModifiedRodriguesParametersElements;
                    unifiedStateModelExponentialMapElementsOutput[ loopIndex ] = convertedUnifiedStateModelExponentialMapElements;

                    // Next step
                    ++loopIndex;
                }
            }
        }
    }

    // Timing
    tend = time( 0 );
    std::cout << "Execution time: " << difftime( tend, tstart ) << " s" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       SAVE OUTPUT                        //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( keplerianElementsInput,
                                          "kepler_input.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< int >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( keplerianElementsUSM7Output,
                                          "kepler_usm7_output.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< int >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( keplerianElementsUSM6Output,
                                          "kepler_usm6_output.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< int >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( keplerianElementsUSMEMOutput,
                                          "kepler_usmem_output.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< int >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( unifiedStateModelQuaternionsElementsOutput,
                                          "usm7_output.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< int >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( unifiedStateModelModifiedRodriguesParametersElementsOutput,
                                          "usm6_output.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< int >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( unifiedStateModelExponentialMapElementsOutput,
                                          "usmem_output.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< int >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
