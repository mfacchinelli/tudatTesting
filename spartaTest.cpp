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
#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"

//! Execute propagation of orbit of Satellite around the Earth.
int main( )
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output;
    using namespace tudat::aerodynamics;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            TABULATED ATMOSPHERE          ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create a tabulated atmosphere object.
    std::map< int, std::string > tabulatedAtmosphereFile = { { 0, getAtmosphereTablesPath( ) + "MCDMeanAtmosphere.dat" } };
    std::vector< AtmosphereDependentVariables > dependentVariables = {
        density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere,
        gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            INPUTS                        ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define inputs
    const std::string SPARTAExecutable = "'/Users/Michele/AE Software/SPARTA/src/spa_mac_mpi'";
    const int numberOfCores = 2;
    const std::string geometryFileUser_ = getSPARTADataPath( ) + "data.sphere"; // check that it is not called data.shape
    const TabulatedAtmosphere atmosphereModel( tabulatedAtmosphereFile, dependentVariables  );
    const std::vector< double > simulationAltitudes = { 100e3 };//{ 100e3, 125e3, 150e3, 200e3, 300e3, 500e3 };
    const int referenceAxis_ = + 0; // axis opposite to free stream velocity and where reference aerodynamic area is taken
    const Eigen::Vector3d momentReferencePoint = Eigen::Vector3d::Zero( );
    const std::string simulationGases_ = "CO2 H O"; // check that it matches gases in SPARTA atmosphere
    const double wallTemperature_ = 300;
    const double accomodationCoefficient_ = 1.0;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            INTERNAL VARIABLES            ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define internal files and paths
    std::string outputDirectory_ = "results";
    std::string outputPath_ = getSPARTADataPath( ) + outputDirectory_;
    std::string inputFile_ = getSPARTADataPath( )  + "in.sparta";
    std::string inputFileTemplate_ = getSPARTADataPath( )  + "SPARTAInputTemplate.txt";
    std::string geometryFileInternal_ = getSPARTADataPath( ) + "data.shape";
    std::string atmosphereSpeciesFile = getSPARTADataPath( ) + "atmosphere.species";
    std::string atmosphereCollisionModelFile = getSPARTADataPath( ) + "atmosphere.vss";

    // Define simulation variables
    double gridSpacing_ = 0.25;
    double simulatedParticlesPerCell_ = 15;
    std::vector< double > simulationAnglesOfAttack = { -75.0, -60.0, -45.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0,
                                                       0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 45.0, 60.0, 75.0 };
    std::vector< double > simulationAnglesOfSideslip = { 0.0 };
    std::vector< double > simulationMolecularSpeedRatios = { 1.0 };
//    { 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0 };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            READ SHAPE FILE               ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Extract simulation box data from shape file

    // Initialize output vectors
    Eigen::Matrix< double, Eigen::Dynamic, 3 > shapePoints_;
    Eigen::Matrix< int, Eigen::Dynamic, 3 > shapeTriangles_;
    int numberOfPoints_ = 0;
    int numberOfTriangles_ = 0;
    {
        // Open file and create file stream.
        std::fstream stream( geometryFileUser_.c_str( ), std::ios::in );

        // Check if file opened correctly.
        if ( stream.fail( ) )
        {
            throw std::runtime_error( "Data file could not be opened: " + geometryFileUser_ );
        }

        // Initialize booleans that specifies once parts of file have been passed.
        bool isNumberOfPointsPassed = false;
        bool isNumberOfTrianglesPassed = false;
        bool isListOfPointsPassed = false;

        // Line based parsing
        std::string line;
        std::vector< std::string > vectorOfIndividualStrings;

        // Read file line-by-line
        int numberOfPointsParsed = 0;
        int numberOfTrianglesParsed = 0;
        while ( !stream.fail( ) && !stream.eof( ) )
        {
            // Get line from stream
            std::getline( stream, line );

            // Trim input string (removes all leading and trailing whitespaces).
            boost::algorithm::trim( line );

            // Skip empty and comment lines
            if ( line.size( ) > 0 && !( line.at( 0 ) == '#' ) )
            {
                // Split string into multiple strings, each containing one element from a line from the data file.
                boost::algorithm::split( vectorOfIndividualStrings,
                                         line,
                                         boost::algorithm::is_any_of( "\t ;, " ),
                                         boost::algorithm::token_compress_on );

                // If this is the first line that is read, it should contain the number of points
                if ( !isNumberOfPointsPassed )
                {
                    if ( vectorOfIndividualStrings.size( ) != 2 )
                    {
                        throw std::runtime_error( "Error when reading multi-array, expected number of points." );
                    }
                    numberOfPoints_ = std::stoi( vectorOfIndividualStrings.at( 0 ) );
                    isNumberOfPointsPassed = true;
                }
                // If this is the first line that is read, it should contain the number of triangles
                else if ( !isNumberOfTrianglesPassed )
                {
                    if ( vectorOfIndividualStrings.size( ) != 2 )
                    {
                        throw std::runtime_error( "Error when reading multi-array, expected number of triangles." );
                    }
                    numberOfTriangles_ = std::stoi( vectorOfIndividualStrings.at( 0 ) );
                    isNumberOfTrianglesPassed = true;
                }
                else if ( !isListOfPointsPassed )
                {
                    if ( vectorOfIndividualStrings.at( 0 ) != "Points" )
                    {
                        // Check line consistency
                        if ( vectorOfIndividualStrings.size( ) != static_cast< unsigned int >( 4 ) )
                        {
                            throw std::runtime_error(
                                        "Error on data line " + std::to_string( numberOfPointsParsed ) +
                                        " found " + std::to_string( vectorOfIndividualStrings.size( ) ) +
                                        " columns, but expected " + std::to_string( 4 ) );
                        }
                        else
                        {
                            // Parse data from current line into output matrix.
                            for ( unsigned int i = 0; i < ( vectorOfIndividualStrings.size( ) - 1 ); i++ )
                            {
                                shapePoints_( numberOfPointsParsed, i ) = std::stod( vectorOfIndividualStrings.at( i + 1 ) );
                            }
                            numberOfPointsParsed++;
                        }
                    }

                    if ( numberOfPointsParsed == numberOfPoints_ )
                    {
                        isListOfPointsPassed = true;
                    }
                }
                else if ( isListOfPointsPassed )
                {
                    if ( vectorOfIndividualStrings.at( 0 ) != "Triangles" )
                    {
                        // Check line consistency
                        if ( vectorOfIndividualStrings.size( ) != static_cast< unsigned int >( 4 ) )
                        {
                            throw std::runtime_error(
                                        "Error on data line " + std::to_string( numberOfTrianglesParsed ) +
                                        " found " + std::to_string( vectorOfIndividualStrings.size( ) ) +
                                        " columns, but expected " + std::to_string( 4 ) );
                        }
                        else
                        {
                            // Parse data from current line into output matrix.
                            for ( unsigned int i = 0; i < ( vectorOfIndividualStrings.size( ) - 1 ); i++ )
                            {
                                shapeTriangles_( numberOfTrianglesParsed, i ) = std::stod( vectorOfIndividualStrings.at( i + 1 ) );
                            }
                            numberOfTrianglesParsed++;
                        }
                    }

                    if ( numberOfTrianglesParsed > numberOfTriangles_ )
                    {
                        throw std::runtime_error( "Number of triangles in file does not match file header." );
                    }
                }

                // Allocate memory for point and triangle matrices.
                if ( isNumberOfPointsPassed && isNumberOfTrianglesPassed )
                {
                    // Check input consistency
                    if ( ( numberOfPoints_ == 0 ) || ( numberOfTriangles_ == 0 ) )
                    {
                        throw std::runtime_error( "Error when reading shape file, expected to find a non-zero number of points and triangles." );
                    }
                    else
                    {
                        // Define size of output vectors
                        shapePoints_.resize( numberOfPoints_, 3 );
                        shapeTriangles_.resize( numberOfTriangles_, 3 );
                    }
                }
            }
        }
    }

    // Get maximum and minimum values in each dimension
    Eigen::Vector3d maximumDimensions_ = shapePoints_.colwise( ).maxCoeff( );
    Eigen::Vector3d minimumDimensions_ = shapePoints_.colwise( ).minCoeff( );
    maximumDimensions_ += 0.5 * maximumDimensions_; // add extra space around shape
    minimumDimensions_ += 0.5 * minimumDimensions_; // add extra space around shape

    // Compute normal to surface elements, area of surface elements and moment arm values
    Eigen::Matrix< double, 3, Eigen::Dynamic > elementSurfaceNormal;
    Eigen::Matrix< double, 1, Eigen::Dynamic > elementSurfaceArea;
    Eigen::Matrix< double, 3, Eigen::Dynamic > elementMomentArm;
    elementSurfaceNormal.resize( 3, numberOfTriangles_ );
    elementSurfaceArea.resize( 1, numberOfTriangles_ );
    elementMomentArm.resize( 3, numberOfTriangles_ );
    Eigen::Matrix3d currentVertices;
    Eigen::Vector3d currentNormal;
    Eigen::Vector3d currentCentroid;
    double currentNormalNorm;
    for ( int i = 0; i < numberOfTriangles_; i++ )
    {
        // Compute properties of current surface element
        for ( unsigned int j = 0; j < 3; j++ )
        {
            currentVertices.row( j ) = shapePoints_.row( shapeTriangles_( i, j ) - 1 );
        }
        currentNormal = ( currentVertices.row( 1 ) - currentVertices.row( 0 ) ).cross(
                    currentVertices.row( 2 ) - currentVertices.row( 0 ) );
        currentNormalNorm = currentNormal.norm( );
        currentCentroid = currentVertices.colwise( ).sum( ) / 3.0;

        // Find normal, area and distance to reference point
        elementSurfaceNormal.col( i ) = currentNormal / currentNormalNorm;
        elementSurfaceArea( i ) = 0.5 * currentNormalNorm;
        elementMomentArm.col( i ) = currentCentroid - momentReferencePoint;
    }

    // Compute cross-sectional area
    Eigen::Vector3d shapeCrossSectionalArea;
    for ( unsigned int i = 0; i < 3; i++ )
    {
        shapeCrossSectionalArea( i ) =
                0.5 * elementSurfaceNormal.row( i ).cwiseAbs( ).dot( elementSurfaceArea );
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            EXTRACT ATMOSPHERE CONDITIONS ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Extract atmosphere conditions (turn into enum or struct?)
    //      0 - density
    //      1 - pressure
    //      2 - temperature
    //      3 - gas constant
    //      4 - molar mass
    //      5 - number density
    std::vector< std::vector< double > > atmosphericConditions_;
    atmosphericConditions_.resize( 6 );
    for ( unsigned int i = 0; i < simulationAltitudes.size( ); i++ )
    {
        atmosphericConditions_.at( 0 ).push_back( atmosphereModel.getDensity( simulationAltitudes.at( i ) ) );
        atmosphericConditions_.at( 1 ).push_back( atmosphereModel.getPressure( simulationAltitudes.at( i ) ) );
        atmosphericConditions_.at( 2 ).push_back( atmosphereModel.getTemperature( simulationAltitudes.at( i ) ) );
        atmosphericConditions_.at( 3 ).push_back( atmosphereModel.getSpecificGasConstant( simulationAltitudes.at( i ) ) );
        atmosphericConditions_.at( 4 ).push_back( tudat::physical_constants::MOLAR_GAS_CONSTANT /
                                                 atmosphericConditions_.at( 3 ).at( i ) );
        atmosphericConditions_.at( 5 ).push_back( tudat::physical_constants::AVOGADRO_CONSTANT *
                                                 atmosphericConditions_.at( 0 ).at( i ) / atmosphericConditions_.at( 4 ).at( i ) );
    }
    // add gases extraction

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            GET SIMULATION CONDITIONS     ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Simulation boundary and grid
    Eigen::Vector6d simulationBoundaries_;
    for ( unsigned int i = 0; i < 3; i++ )
    {
        simulationBoundaries_( 2 * i ) = minimumDimensions_( i );
        simulationBoundaries_( 2 * i + 1 ) = maximumDimensions_( i );
    }
    Eigen::Vector3d simulationGrid_ = ( maximumDimensions_ - minimumDimensions_ ) / gridSpacing_;

    // Convert molecular speed ratio to stream velocity and compute simulation time step and ratio of real to simulated variables
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > freeStreamVelocities_;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > simulationTimeStep_;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > ratioOfRealToSimulatedParticles_;
    freeStreamVelocities_.resize( simulationAltitudes.size( ), simulationMolecularSpeedRatios.size( ) );
    simulationTimeStep_.resize( simulationAltitudes.size( ), simulationMolecularSpeedRatios.size( ) );
    ratioOfRealToSimulatedParticles_.resize( simulationAltitudes.size( ), 1 );
    for ( unsigned int h = 0; h < simulationAltitudes.size( ); h++ )
    {
        for ( unsigned int s = 0; s < simulationMolecularSpeedRatios.size( ); s++ )
        {
            freeStreamVelocities_( h, s ) = simulationMolecularSpeedRatios.at( s ) * std::sqrt(
                        2.0 * atmosphericConditions_.at( 3 ).at( h ) * atmosphericConditions_.at( 2 ).at( h ) );
            simulationTimeStep_( h, s ) = 0.1 * ( maximumDimensions_( std::abs( referenceAxis_ ) ) -
                                                 minimumDimensions_( std::abs( referenceAxis_ ) ) ) /
                    freeStreamVelocities_( h, s );
            // time step is taken as time it takes for a particle to travel for 10 % of the box
        }
        ratioOfRealToSimulatedParticles_( h ) = atmosphericConditions_.at( 5 ).at( h ) *
                std::pow( gridSpacing_, 3 ) / simulatedParticlesPerCell_;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            PROCESS FILES                 ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Initialize output variable
    std::string inputTemplate_;

    // Open file and create file stream.
    {
        std::fstream stream( inputFileTemplate_.c_str( ), std::ios::in );

        // Check if file opened correctly.
        if ( stream.fail( ) )
        {
            throw std::runtime_error( "Data file could not be opened: " + inputFileTemplate_ );
        }

        // Line based parsing
        std::string line;

        // Read file line-by-line
        while ( !stream.fail( ) && !stream.eof( ) )
        {
            // Get line from stream
            std::getline( stream, line );

            // Trim input string (removes all leading and trailing whitespaces).
            boost::algorithm::trim( line );

            // Skip empty and comment lines
            if ( line.size( ) > 0 && !( line.at( 0 ) == '#' ) )
            {
                // Append line to string
                inputTemplate.append( line + "\n" );
            }
        }
    }

    // Copy input shape file to default name
    std::string commandString = "cp " + geometryFileUser_ + " " + geometryFileInternal_;
    std::system( commandString.c_str( ) );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            START SPARTA ANALYSIS         ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Generate command string for SPARTA
    std::cout << "Initiating SPARTA simulation. This may take a while." << std::endl;
    std::string runSPARTACommandString = "cd " + getSPARTADataPath( ) + "; echo off; " +
            SPARTAExecutable +
            " -in " + inputFile_ + "; echo on";
    // "mpirun -np " + std::to_string( numberOfCores ) + " " +

    // Loop over simulation parameters and run SPARTA
    boost::multi_array< Eigen::Vector6d, 3 > aerodynamicCoefficients_( boost::extents[ simulationAltitudes.size( ) ][
                simulationMolecularSpeedRatios.size( ) ][ simulationAnglesOfAttack.size( ) ] );
    std::string anglesOfAttack;
    Eigen::Vector3d velocityVector;
    for ( unsigned int h = 0; h < simulationAltitudes.size( ); h++ )
    {
        for ( unsigned int s = 0; s < simulationMolecularSpeedRatios.size( ); s++ )
        {
            // Get velocity vector
            velocityVector = Eigen::Vector3d::Zero( );
            velocityVector( std::abs( referenceAxis_ ) ) = ( std::signbit( referenceAxis_ ) ? 1.0 : -1.0 ) *
                    freeStreamVelocities_( h, s );

            // Get angles of attack string
            for ( double a : simulationAnglesOfAttack )
            {
                anglesOfAttack += printToStringWithPrecision( a, 0 ) + " ";
            }

            // Print to file
            FILE * fileIdentifier = std::fopen( inputFile_.c_str( ), "w" );
            std::fprintf( fileIdentifier, inputTemplate_.c_str( ), simulationBoundaries_( 0 ), simulationBoundaries_( 1 ),
                          simulationBoundaries_( 2 ), simulationBoundaries_( 3 ), simulationBoundaries_( 4 ),
                          simulationBoundaries_( 5 ), simulationGrid_( 0 ), simulationGrid_( 1 ), simulationGrid_( 2 ),
                          atmosphericConditions_.at( 5 ).at( h ), ratioOfRealToSimulatedParticles_( h ), simulationGases_.c_str( ),
                          simulationGases_.c_str( ), velocityVector( 0 ), velocityVector( 1 ), velocityVector( 2 ),
                          simulationGases_.c_str( ), atmosphericConditions_.at( 2 ).at( h ), anglesOfAttack.c_str( ),
                          wallTemperature_, accomodationCoefficient_, simulationTimeStep_( h, s ), outputDirectory_.c_str( ) );
            std::fclose( fileIdentifier );

            // Run SPARTA
            int systemStatus = std::system( runSPARTACommandString.c_str( ) );
            if ( systemStatus != 0 )
            {
                throw std::runtime_error( "Error: SPARTA simulation failed. See the log.sparta file for more details." );
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////            PROCESS RESULTS               ////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Loop over angles of attack
            std::string temporaryOutputFile;
            std::vector< std::string > outputFileExtensions = { ".400", ".600", ".800", ".1000" };
            Eigen::Matrix< double, Eigen::Dynamic, 7 > outputMatrix;
            Eigen::Matrix< double, 3, Eigen::Dynamic > meanPressureValues;
            Eigen::Matrix< double, 3, Eigen::Dynamic > meanShearValues;
            meanPressureValues.resize( 3, numberOfTriangles_ );
            meanShearValues.resize( 3, numberOfTriangles_ );
            for ( unsigned int a = 0; a < simulationAnglesOfAttack.size( ); a++ )
            {
                // Get file name
                temporaryOutputFile = outputPath_ + "/" + printToStringWithPrecision(
                            simulationAnglesOfAttack.at( a ), 0 ) + ".coeff";

                // Read output files and compute mean pressure and shear force values
                meanPressureValues.setZero( );
                meanShearValues.setZero( );
                for ( unsigned int i = 0; i < outputFileExtensions.size( ); i++ )
                {
                    outputMatrix = readMatrixFromFile( temporaryOutputFile + outputFileExtensions.at( i ), "\t ;,", "%", 9 );
                    for ( unsigned int j = 0; j < 3; j++ )
                    {
                        meanPressureValues.row( j ) += outputMatrix.col( j + 1 ).transpose( );
                        meanShearValues.row( j ) += outputMatrix.col( j + 4 ).transpose( );
                    }
               }
                meanPressureValues /= outputFileExtensions.size( );
                meanShearValues /= outputFileExtensions.size( );

                // Convert pressure and shear forces to coefficients
                aerodynamicCoefficients_[ h ][ s ][ a ] = computeAerodynamicCoefficientsFromPressureShear(
                            meanPressureValues,
                            meanShearValues,
                            atmosphericConditions_.at( 0 ).at( h ), // density
                            freeStreamVelocities_( h, s ),
                            atmosphericConditions_.at( 1 ).at( h ), // pressure
                            elementSurfaceNormal,
                            elementSurfaceArea,
                            elementMomentArm,
                            shapeCrossSectionalArea( std::abs( referenceAxis_ ) ) );
                std::cout << std::endl << "Altitude: " << simulationAltitudes.at( h ) << std::endl
                          << "Speed Ratio: " << simulationMolecularSpeedRatios.at( s ) << std::endl
                          << "Angle of Attack: " << simulationAnglesOfAttack.at( a ) << std::endl
                          << "Coefficients: " << aerodynamicCoefficients_[ h ][ s ][ a ].transpose( ) << std::endl;
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Here." << std::endl;

    // Clean up results folder
    commandString = "rm " + outputPath_ + "/*"; // overwrite
    std::system( commandString.c_str( ) );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
