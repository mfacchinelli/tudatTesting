/*    Copyright (c) 2010-2018, Delft University of Technology
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

    // Select mode
    bool fullAngleSpectrum = true;

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
            inclinationAngles.push_back( convertDegreesToRadians( angle ) );
        }
        else if ( angle == 180.0 )
        {
            inclinationAngles.push_back( convertDegreesToRadians( angle ) - 0.001 );
        }

        // Add other angles up to 360 degrees
        rightAscensionOfAscendingNodeAngles.push_back( convertDegreesToRadians( angle ) );
        argumentOfPeriapsisAngles.push_back( convertDegreesToRadians( angle ) );
        trueAnomalyAngles.push_back( convertDegreesToRadians( angle ) );

        // Next step
        angle += fullAngleSpectrum ? 15.0 : 360.0;
    }

    // Create Keplerian elements vector
    Eigen::Vector6d keplerianElements = Eigen::Vector6d::Zero( );
    keplerianElements( semiMajorAxisIndex ) = 6378.1363e3 + 1000e3;
    keplerianElements( eccentricityIndex ) = 0.25;

    // Standard gravitational parameter
    const double centralBodyGravitationalParameter = 3.9860044189e14; // [m^3/s^2]

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CONVERT TO AND FROM USM            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create map of results
    std::map< int, Eigen::Vector6d > keplerianElementsInput;
    std::map< int, Eigen::Vector6d > keplerianElementsUSM7Output;
    std::map< int, Eigen::Vector6d > keplerianElementsUSM6Output;
    std::map< int, Eigen::Vector6d > keplerianElementsUSMEMOutput;
    std::map< int, Eigen::Vector7d > unifiedStateModelQuaternionsElementsOutput;
    std::map< int, Eigen::Vector7d > unifiedStateModelModifiedRodriguesParametersElementsOutput;
    std::map< int, Eigen::Vector7d > unifiedStateModelExponentialMapElementsOutput;

    // Predefine Keplerian vectors
    Eigen::Vector6d convertedKeplerianElementsUSM7;
    Eigen::Vector6d convertedKeplerianElementsUSM6;
    Eigen::Vector6d convertedKeplerianElementsUSMEM;

    // Predefine Cartesian vectors
    Eigen::Vector6d cartesianElements;
    Eigen::Vector6d convertedCartesianElementsUSM7;
    Eigen::Vector6d convertedCartesianElementsUSM6;
    Eigen::Vector6d convertedCartesianElementsUSMEM;

    // Predefine USM vectors
    Eigen::Vector7d convertedUnifiedStateModelQuaternionsElements;
    Eigen::Vector7d convertedUnifiedStateModelModifiedRodriguesParameterElements;
    Eigen::Vector7d convertedUnifiedStateModelExponentialMapElements;

    // Select conversion mode
    //      0: Kepl -> USM -> Kepl
    //      1: Kepl -> Cart -> USM -> Cart -> Kepl
    //      2: Kepl -> Cart -> USM -> Kepl
    //      3: Kepl -> USM -> Cart -> Kepl
    bool fullConversionAnalysis = true;
    int conversionInitialValue = fullConversionAnalysis ? 0 : 0;
    int conversionFinalValue = fullConversionAnalysis ? 4 : ( conversionInitialValue + 1 );

    // Start loop
    for ( int conversionMode = conversionInitialValue; conversionMode < conversionFinalValue; conversionMode++ )
    {
        // Print mode
        std::cout << "Mode: " << conversionMode << std::endl;

        // Timeing
        time_t tstart, tend;
        tstart = time( 0 );

        // Loop over angles
        int loopIndex = 1;
        double tolerance = 1e-15;
        for ( double i : inclinationAngles )
        {
            for ( double O : rightAscensionOfAscendingNodeAngles )
            {
                for ( double o : argumentOfPeriapsisAngles )
                {
                    for ( double t : trueAnomalyAngles )
                    {
                        // Complete Keplerian elements vector
                        if ( fullAngleSpectrum )
                        {
                            keplerianElements( inclinationIndex ) = i;
                            keplerianElements( longitudeOfAscendingNodeIndex ) = O;
                            keplerianElements( argumentOfPeriapsisIndex ) = o;
                            keplerianElements( trueAnomalyIndex ) = t;
                        }
                        else
                        {
                            keplerianElements( inclinationIndex ) = convertDegreesToRadians< double >( 179.9427042204868 );
                            keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians< double >( 0 );
                            keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians< double >( 30 );
                            keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians< double >( 60 );
                            std::cout << "Kepler In: " << keplerianElements.transpose( ) << std::endl;
                        }

                        // Check for inconsistencies
                        if ( std::abs( keplerianElements( inclinationIndex ) ) < tolerance )
                        {
                            keplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
                        }

                        // Convert directly to/from Keplerian or pass through Cartesian, as well
                        switch ( conversionMode )
                        {
                        case 0:
                        {
                            // Convert to USM7 and back
                            convertedUnifiedStateModelQuaternionsElements =
                                    convertKeplerianToUnifiedStateModelQuaternionsElements(
                                        keplerianElements, centralBodyGravitationalParameter );

                            convertedKeplerianElementsUSM7 = convertUnifiedStateModelQuaternionsToKeplerianElements(
                                        convertedUnifiedStateModelQuaternionsElements, centralBodyGravitationalParameter );

                            // Convert to USM6 and back
                            convertedUnifiedStateModelModifiedRodriguesParameterElements =
                                    convertKeplerianToUnifiedStateModelModifiedRodriguesParameterElements(
                                        keplerianElements, centralBodyGravitationalParameter );

                            convertedKeplerianElementsUSM6 =
                                    convertUnifiedStateModelModifiedRodriguesParametersToKeplerianElements(
                                        convertedUnifiedStateModelModifiedRodriguesParameterElements,
                                        centralBodyGravitationalParameter );

                            // Convert to USMEM and back
                            convertedUnifiedStateModelExponentialMapElements =
                                    convertKeplerianToUnifiedStateModelExponentialMapElements(
                                        keplerianElements, centralBodyGravitationalParameter );

                            convertedKeplerianElementsUSMEM = convertUnifiedStateModelExponentialMapToKeplerianElements(
                                        convertedUnifiedStateModelExponentialMapElements, centralBodyGravitationalParameter );
                            break;
                        }
                        case 1:
                        {
                            // Get Cartesian elements
                            cartesianElements = convertKeplerianToCartesianElements(
                                        keplerianElements, centralBodyGravitationalParameter );

                            // Convert to USM7 and back
                            convertedUnifiedStateModelQuaternionsElements = convertCartesianToUnifiedStateModelQuaternionsElements(
                                        cartesianElements, centralBodyGravitationalParameter );

                            convertedCartesianElementsUSM7 = convertUnifiedStateModelQuaternionsToCartesianElements(
                                        convertedUnifiedStateModelQuaternionsElements, centralBodyGravitationalParameter );
                            convertedKeplerianElementsUSM7 = convertCartesianToKeplerianElements(
                                        convertedCartesianElementsUSM7, centralBodyGravitationalParameter );

                            // Convert to USM6 and back
                            convertedUnifiedStateModelModifiedRodriguesParameterElements =
                                    convertCartesianToUnifiedStateModelModifiedRodriguesParameterElements(
                                        cartesianElements, centralBodyGravitationalParameter );

                            convertedCartesianElementsUSM6 =
                                    convertUnifiedStateModelModifiedRodriguesParametersToCartesianElements(
                                        convertedUnifiedStateModelModifiedRodriguesParameterElements,
                                        centralBodyGravitationalParameter );
                            convertedKeplerianElementsUSM6 = convertCartesianToKeplerianElements(
                                        convertedCartesianElementsUSM6, centralBodyGravitationalParameter );

                            // Convert to USMEM and back
                            convertedUnifiedStateModelExponentialMapElements =
                                    convertCartesianToUnifiedStateModelExponentialMapElements(
                                        cartesianElements, centralBodyGravitationalParameter );

                            convertedCartesianElementsUSMEM =
                                    convertUnifiedStateModelExponentialMapToCartesianElements(
                                        convertedUnifiedStateModelExponentialMapElements,
                                        centralBodyGravitationalParameter );
                            convertedKeplerianElementsUSMEM = convertCartesianToKeplerianElements(
                                        convertedCartesianElementsUSMEM, centralBodyGravitationalParameter );
                            break;
                        }
                        case 2:
                        {
                            // Get Cartesian elements
                            cartesianElements = convertKeplerianToCartesianElements(
                                        keplerianElements, centralBodyGravitationalParameter );

                            // Convert to USM7 and back
                            convertedUnifiedStateModelQuaternionsElements = convertCartesianToUnifiedStateModelQuaternionsElements(
                                        cartesianElements, centralBodyGravitationalParameter );

                            convertedKeplerianElementsUSM7 = convertUnifiedStateModelQuaternionsToKeplerianElements(
                                        convertedUnifiedStateModelQuaternionsElements, centralBodyGravitationalParameter );

                            // Convert to USM6 and back
                            convertedUnifiedStateModelModifiedRodriguesParameterElements =
                                    convertCartesianToUnifiedStateModelModifiedRodriguesParameterElements(
                                        cartesianElements, centralBodyGravitationalParameter );

                            convertedKeplerianElementsUSM6 =
                                    convertUnifiedStateModelModifiedRodriguesParametersToKeplerianElements(
                                        convertedUnifiedStateModelModifiedRodriguesParameterElements,
                                        centralBodyGravitationalParameter );

                            // Convert to USMEM and back
                            convertedUnifiedStateModelExponentialMapElements =
                                    convertCartesianToUnifiedStateModelExponentialMapElements(
                                        cartesianElements, centralBodyGravitationalParameter );

                            convertedKeplerianElementsUSMEM =
                                    convertUnifiedStateModelExponentialMapToKeplerianElements(
                                        convertedUnifiedStateModelExponentialMapElements,
                                        centralBodyGravitationalParameter );
                            break;
                        }
                        case 3:
                        {
                            // Convert to USM7 and back
                            convertedUnifiedStateModelQuaternionsElements = convertKeplerianToUnifiedStateModelQuaternionsElements(
                                        keplerianElements, centralBodyGravitationalParameter );

                            convertedCartesianElementsUSM7 = convertUnifiedStateModelQuaternionsToCartesianElements(
                                        convertedUnifiedStateModelQuaternionsElements, centralBodyGravitationalParameter );
                            convertedKeplerianElementsUSM7 = convertCartesianToKeplerianElements(
                                        convertedCartesianElementsUSM7, centralBodyGravitationalParameter );

                            // Convert to USM6 and back
                            convertedUnifiedStateModelModifiedRodriguesParameterElements =
                                    convertKeplerianToUnifiedStateModelModifiedRodriguesParameterElements(
                                        keplerianElements, centralBodyGravitationalParameter );

                            convertedCartesianElementsUSM6 =
                                    convertUnifiedStateModelModifiedRodriguesParametersToCartesianElements(
                                        convertedUnifiedStateModelModifiedRodriguesParameterElements,
                                        centralBodyGravitationalParameter );
                            convertedKeplerianElementsUSM6 = convertCartesianToKeplerianElements(
                                        convertedCartesianElementsUSM6, centralBodyGravitationalParameter );

                            // Convert to USMEM and back
                            convertedUnifiedStateModelExponentialMapElements =
                                    convertKeplerianToUnifiedStateModelExponentialMapElements(
                                        keplerianElements, centralBodyGravitationalParameter );

                            convertedCartesianElementsUSMEM =
                                    convertUnifiedStateModelExponentialMapToCartesianElements(
                                        convertedUnifiedStateModelExponentialMapElements,
                                        centralBodyGravitationalParameter );
                            convertedKeplerianElementsUSMEM = convertCartesianToKeplerianElements(
                                        convertedCartesianElementsUSMEM, centralBodyGravitationalParameter );
                            break;
                        }
                        default:
                            throw std::runtime_error( "Conversion mode not recognized." );
                        }

                        // Save to map
                        keplerianElementsInput[ loopIndex ] = keplerianElements;
                        keplerianElementsUSM7Output[ loopIndex ] = convertedKeplerianElementsUSM7;
                        keplerianElementsUSM6Output[ loopIndex ] = convertedKeplerianElementsUSM6;
                        keplerianElementsUSMEMOutput[ loopIndex ] = convertedKeplerianElementsUSMEM;
                        unifiedStateModelQuaternionsElementsOutput[ loopIndex ] = convertedUnifiedStateModelQuaternionsElements;
                        unifiedStateModelModifiedRodriguesParametersElementsOutput[ loopIndex ] = convertedUnifiedStateModelModifiedRodriguesParameterElements;
                        unifiedStateModelExponentialMapElementsOutput[ loopIndex ] = convertedUnifiedStateModelExponentialMapElements;

                        // Next step
                        ++loopIndex;
                    }
                }
            }
        }

        // Timing
        tend = time( 0 );
        std::cout << "Conversion time: " << difftime( tend, tstart ) << " s" << std::endl;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////       SAVE OUTPUT                        //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Write perturbed satellite propagation history to file.
        writeDataMapToTextFile( keplerianElementsInput,
                                "kepler_input_" + std::to_string( conversionMode ) + ".dat", getOutputPath( "Conversions" ) );

        writeDataMapToTextFile( keplerianElementsUSM7Output,
                                "kepler_usm7_output_" + std::to_string( conversionMode ) + ".dat", getOutputPath( "Conversions" ) );

        writeDataMapToTextFile( keplerianElementsUSM6Output,
                                "kepler_usm6_output_" + std::to_string( conversionMode ) + ".dat", getOutputPath( "Conversions" ) );

        writeDataMapToTextFile( keplerianElementsUSMEMOutput,
                                "kepler_usmem_output_" + std::to_string( conversionMode ) + ".dat", getOutputPath( "Conversions" ) );

        writeDataMapToTextFile( unifiedStateModelQuaternionsElementsOutput,
                                "usm7_output_" + std::to_string( conversionMode ) + ".dat", getOutputPath( "Conversions" ) );

        writeDataMapToTextFile( unifiedStateModelModifiedRodriguesParametersElementsOutput,
                                "usm6_output_" + std::to_string( conversionMode ) + ".dat", getOutputPath( "Conversions" ) );

        writeDataMapToTextFile( unifiedStateModelExponentialMapElementsOutput,
                                "usmem_output_" + std::to_string( conversionMode ) + ".dat", getOutputPath( "Conversions" ) );
    }

    // Final statement.
    return EXIT_SUCCESS;
}
