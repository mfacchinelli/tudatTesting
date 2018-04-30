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

template< typename T >
std::string to_string_with_precision( const T a_value, const int n = 3 )
{
    std::ostringstream out;
    out << std::setprecision( n ) << std::fixed << a_value;
    return out.str( );
}

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
    std::string SPARTAExecutable = "'/Users/Michele/AE Software/SPARTA/src/spa_mac_mpi'";
    int numberOfCores = 2;
    std::string geometryFileUser = getSPARTADataPath( ) + "data.sphere"; // check that it is not called data.shape
    TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFile, dependentVariables  );
    std::vector< double > simulationAltitudes = { 100e3, 125e3, 150e3, 200e3, 300e3, 500e3 };
    int referenceAxis = + 0; // axis opposite to free stream velocity and where reference aerodynamic area is taken
    std::string simulationGases = "CO2 H O"; // check that it matches gases in SPARTA atmosphere
    double atmosphereTemperature = 200;
    double wallTemperature = 300;
    double accomodationCoefficient = 1.0;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            INTERNAL VARIABLES            ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define internal files and paths
    std::string outputDirectory = "results";
    std::string outputPath = getSPARTADataPath( ) + outputDirectory;
    std::string inputFile = getSPARTADataPath( )  + "in.sparta";
    std::string inputFileTemplate = getSPARTADataPath( )  + "SPARTAInputTemplate.txt";
    std::string geometryFileInternal = getSPARTADataPath( ) + "data.shape";
    std::string atmosphereSpeciesFile = getSPARTADataPath( ) + "atmosphere.species";
    std::string atmosphereCollisionModelFile = getSPARTADataPath( ) + "atmosphere.vss";

    // Define simulation variables
    double gridSpacing = 0.25;
    double simulatedParticlesPerCell = 15;
    std::vector< double > simulationAnglesOfAttack = { -75.0, -60.0, -45.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0,
                                                       0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 45.0, 60.0, 75.0 };
    std::vector< double > simulationAnglesOfSideslip = { 0.0 };
    std::vector< double > simulationMolecularSpeedRatios = { 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0 };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            READ SHAPE FILE               ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Extract simulation box data from shape file

    // Initialize output vectors
    Eigen::Matrix< double, Eigen::Dynamic, 3 > shapePoints;
    Eigen::Matrix< int, Eigen::Dynamic, 3 > shapeTriangles;
    int numberOfPoints = 0;
    int numberOfTriangles = 0;

    {
        // Open file and create file stream.
        std::fstream stream( geometryFileUser.c_str( ), std::ios::in );

        // Check if file opened correctly.
        if ( stream.fail( ) )
        {
            throw std::runtime_error( "Data file could not be opened: " + geometryFileUser );
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
                    numberOfPoints = std::stoi( vectorOfIndividualStrings.at( 0 ) );
                    isNumberOfPointsPassed = true;
                }
                // If this is the first line that is read, it should contain the number of triangles
                else if ( !isNumberOfTrianglesPassed )
                {
                    if ( vectorOfIndividualStrings.size( ) != 2 )
                    {
                        throw std::runtime_error( "Error when reading multi-array, expected number of triangles." );
                    }
                    numberOfTriangles = std::stoi( vectorOfIndividualStrings.at( 0 ) );
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
                                shapePoints( numberOfPointsParsed, i ) = std::stod( vectorOfIndividualStrings.at( i + 1 ) );
                            }
                            numberOfPointsParsed++;
                        }
                    }

                    if ( numberOfPointsParsed == numberOfPoints )
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
                                shapeTriangles( numberOfTrianglesParsed, i ) = std::stod( vectorOfIndividualStrings.at( i + 1 ) );
                            }
                            numberOfTrianglesParsed++;
                        }
                    }

                    if ( numberOfTrianglesParsed > numberOfTriangles )
                    {
                        throw std::runtime_error( "Number of triangles in file does not match file header." );
                    }
                }

                // Allocate memory for point and triangle matrices.
                if ( isNumberOfPointsPassed && isNumberOfTrianglesPassed )
                {
                    // Check input consistency
                    if ( ( numberOfPoints == 0 ) || ( numberOfTriangles == 0 ) )
                    {
                        throw std::runtime_error( "Error when reading shape file, expected to find a non-zero number of points and triangles." );
                    }
                    else
                    {
                        // Define size of output vectors
                        shapePoints.resize( numberOfPoints, 3 );
                        shapeTriangles.resize( numberOfTriangles, 3 );
                    }
                }
            }
        }
    }

    // Get maximum and minimum values in each dimension
    Eigen::Vector3d maximumDimensions = shapePoints.colwise( ).maxCoeff( );
    Eigen::Vector3d minimumDimensions = shapePoints.colwise( ).minCoeff( );
    maximumDimensions += 0.5 * maximumDimensions; // add extra space around shape
    minimumDimensions += 0.5 * minimumDimensions; // add extra space around shape

    // Compute normal to surface elements, area of surface elements and cross-sectional area
    Eigen::Matrix< double, Eigen::Dynamic, 3 > elementSurfaceNormal;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > elementSurfaceArea;
    Eigen::Vector3d shapeCrossSectionalArea;
    elementSurfaceNormal.resize( numberOfTriangles, 3 );
    elementSurfaceArea.resize( numberOfTriangles, 3 );
    Eigen::Matrix3d currentVertices;
    Eigen::Vector3d currentNormal;
    double currentNormalNorm;
    for ( int i = 0; i < numberOfTriangles; i++ )
    {
        for ( unsigned int j = 0; j < 3; j++ )
        {
            currentVertices.row( j ) = shapePoints.row( shapeTriangles( i, j ) - 1 );
        }
        currentNormal = ( currentVertices.row( 1 ) - currentVertices.row( 0 ) ).cross(
                    currentVertices.row( 2 ) - currentVertices.row( 0 ) );
        currentNormalNorm = currentNormal.norm( );
        elementSurfaceArea( i ) = 0.5 * currentNormalNorm;
        elementSurfaceNormal.row( i ) = currentNormal / currentNormalNorm;
    }
    for ( unsigned int i = 0; i < 3; i++ )
    {
        shapeCrossSectionalArea( i ) =
                0.5 * elementSurfaceNormal.col( i ).cwiseAbs( ).transpose( ).dot( elementSurfaceArea );
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            EXTRACT ATMOSPHERE CONDITIONS ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Extract atmosphere conditions (turn into enum)
    //      0 - density
    //      1 - temperature
    //      2 - gas constant
    //      3 - molar mass
    //      4 - number density
    std::vector< std::vector< double > > atmosphericConditions;
    atmosphericConditions.resize( 5 );
    for ( unsigned int i = 0; i < simulationAltitudes.size( ); i++ )
    {
        atmosphericConditions.at( 0 ).push_back( tabulatedAtmosphere.getDensity( simulationAltitudes.at( i ) ) );
        atmosphericConditions.at( 1 ).push_back( tabulatedAtmosphere.getTemperature( simulationAltitudes.at( i ) ) );
        atmosphericConditions.at( 2 ).push_back( tabulatedAtmosphere.getSpecificGasConstant( simulationAltitudes.at( i ) ) );
        atmosphericConditions.at( 3 ).push_back( tudat::physical_constants::MOLAR_GAS_CONSTANT /
                                                 atmosphericConditions.at( 2 ).at( i ) );
        atmosphericConditions.at( 4 ).push_back( tudat::physical_constants::AVOGADRO_CONSTANT *
                                                 atmosphericConditions.at( 0 ).at( i ) / atmosphericConditions.at( 3 ).at( i ) );
    }
    // add gases extraction

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            GET SIMULATION CONDITIONS     ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Simulation boundary and grid
    Eigen::Vector6d simulationBoundaries;
    for ( unsigned int i = 0; i < 3; i++ )
    {
        simulationBoundaries( 2 * i ) = minimumDimensions( i );
        simulationBoundaries( 2 * i + 1 ) = maximumDimensions( i );
    }
    Eigen::Vector3d simulationGrid = ( maximumDimensions - minimumDimensions ) / gridSpacing;

    // Convert molecular speed ratio to stream velocity and compute simulation time step and ratio of real to simulated variables
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > freeStreamVelocities;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > simulationTimeStep;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > ratioOfRealToSimulatedParticles;
    freeStreamVelocities.resize( simulationAltitudes.size( ), simulationMolecularSpeedRatios.size( ) );
    simulationTimeStep.resize( simulationAltitudes.size( ), simulationMolecularSpeedRatios.size( ) );
    ratioOfRealToSimulatedParticles.resize( simulationAltitudes.size( ), 1 );
    for ( unsigned int h = 0; h < simulationAltitudes.size( ); h++ )
    {
        for ( unsigned int s = 0; s < simulationMolecularSpeedRatios.size( ); s++ )
        {
            freeStreamVelocities( h, s ) = simulationMolecularSpeedRatios.at( s ) * std::sqrt(
                        2.0 * atmosphericConditions.at( 2 ).at( h ) * atmosphericConditions.at( 1 ).at( h ) );
            simulationTimeStep( h, s ) = 0.1 * ( maximumDimensions( std::abs( referenceAxis ) ) -
                                                 minimumDimensions( std::abs( referenceAxis ) ) ) /
                    freeStreamVelocities( h, s );
            // time step is taken as time it takes for a particle to travel for 10 % of the box
        }
        ratioOfRealToSimulatedParticles( h ) = atmosphericConditions.at( 4 ).at( h ) *
                std::pow( gridSpacing, 3 ) / simulatedParticlesPerCell;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            PROCESS FILES                 ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Initialize output variable
    std::string inputTemplate;

    // Open file and create file stream.
    {
        std::fstream stream( inputFileTemplate.c_str( ), std::ios::in );

        // Check if file opened correctly.
        if ( stream.fail( ) )
        {
            throw std::runtime_error( "Data file could not be opened: " + inputFileTemplate );
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
    std::string commandString = "cp " + geometryFileUser + " " + geometryFileInternal;
    std::system( commandString.c_str( ) );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            START SPARTA ANALYSIS         ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Warn user and delete all files
    std::cout << "Warning in rarefied flow analysis with SPARTA, running the "
                 "simulation will delete all existing results" << std::endl;
    commandString = "rm " + outputPath + "/*/*.coeff.*"; // overwrite
    std::system( commandString.c_str( ) );

    // Generate command string for SPARTA
    std::string runSPARTACommandString = "cd " + getSPARTADataPath( ) + "; " +
            SPARTAExecutable +
            " -in " + inputFile;
    // "mpirun -np " + std::to_string( numberOfCores ) + " " +

    // Loop over simulation parameters and run SPARTA
    std::string anglesOfAttack;
    Eigen::Vector3d velocityVector;
    std::string temporaryOutputDirectory;
    for ( unsigned int h = 0; h < simulationAltitudes.size( ); h++ )
    {
        for ( unsigned int s = 0; s < simulationMolecularSpeedRatios.size( ); s++ )
        {

        }
    }
    {
        unsigned int h = 0;
        unsigned int s = 0;

        // Get velocity vector
        velocityVector = Eigen::Vector3d::Zero( );
        velocityVector( std::abs( referenceAxis ) ) = ( std::signbit( referenceAxis ) ? 1.0 : -1.0 ) *
                freeStreamVelocities( h, s );

        // Get angles of attack string
        for ( double a : simulationAnglesOfAttack )
        {
            anglesOfAttack += to_string_with_precision< double >( a, 2 ) + " ";
        }

        // Get output directory for these conditions
        temporaryOutputDirectory = outputDirectory + "/" +
                to_string_with_precision< double >( simulationAltitudes.at( h ) / 1e3, 0 ) + "/" +
                to_string_with_precision< double >( simulationMolecularSpeedRatios.at( s ), 0 );
        std::cout << temporaryOutputDirectory << std::endl;

        // Print to file
        FILE * fileIdentifier = std::fopen( inputFile.c_str( ), "w" );
        std::fprintf( fileIdentifier, inputTemplate.c_str( ), simulationBoundaries( 0 ), simulationBoundaries( 1 ),
                      simulationBoundaries( 2 ), simulationBoundaries( 3 ), simulationBoundaries( 4 ),
                      simulationBoundaries( 5 ), simulationGrid( 0 ), simulationGrid( 1 ), simulationGrid( 2 ),
                      atmosphericConditions.at( 4 ).at( h ), ratioOfRealToSimulatedParticles( h ), simulationGases.c_str( ),
                      simulationGases.c_str( ), velocityVector( 0 ), velocityVector( 1 ), velocityVector( 2 ),
                      simulationGases.c_str( ), atmosphereTemperature, anglesOfAttack.c_str( ), wallTemperature,
                      accomodationCoefficient, simulationTimeStep( h, s ), temporaryOutputDirectory.c_str( ) );
        std::fclose( fileIdentifier );

        // Run SPARTA
        int systemStatus = std::system( runSPARTACommandString.c_str( ) );
        if ( systemStatus != 0 )
        {
            throw std::runtime_error( "Error: SPARTA simulation failed. See log.sparta for more details." );
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            PROCESS RESULTS               ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
