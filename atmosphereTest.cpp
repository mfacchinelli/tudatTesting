/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/InputOutput/basicInputOutput.h"

//! Execute propagation of orbit of Asterix around Mars.
int main( )
{
    using namespace tudat;
    using namespace tudat::aerodynamics;
    using namespace tudat::unit_conversions;

    // Declare tolerance used for Boost tests.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Select mode
    //      0: Tabulated atmosphere test
    //      1: Multi-array reader test
    const int mode = 0;
    switch ( mode )
    {
    case 0:
    {
        // Declare input and output
        double altitudeInput, longitudeInput, latitudeInput, densityInput, pressureInput, temperatureInput;
        double densityOutput, pressureOutput, temperatureOutput;

        // Select case
        //      0: 1-D: corner conditions, IND: nominal order, DEP: nominal order
        //      1: 1-D: interpolated conditions, IND: nominal order, DEP: nominal order
        //      2: 1-D: corner conditions, IND: shuffled order, DEP: shuffled order
        //      3: 3-D: corner conditions, IND: shuffled order, DEP: nominal order
        //      4: 3-D: non-interpolated conditions, IND: shuffled order, DEP: shuffled order
        //      5: 3-D: interpolated conditions, IND: shuffled order, DEP: nominal order
        const int testCase = 4;
        std::cout << "Test Case < " << testCase << " > Selected" << std::endl;
        std::cout << "Tolerance is set at: " << tolerance << std::endl << std::endl;
        switch ( testCase )
        {
        case 0: // 1-D: corner conditions, IND: nominal order, DEP: nominal order
        {
            // Create a tabulated atmosphere object.
            std::string tabulatedAtmosphereFile = input_output::getAtmosphereTablesPath( ) +
                    "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
            TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFile );

            // Declare input and output
            altitudeInput = 0.0;
            densityInput = 1.225;
            pressureInput = 101320.0;
            temperatureInput = 288.15;

            densityOutput = tabulatedAtmosphere.getDensity( altitudeInput );
            pressureOutput = tabulatedAtmosphere.getPressure( altitudeInput );
            temperatureOutput = tabulatedAtmosphere.getTemperature( altitudeInput );
            break;
        }
        case 1: // 1-D: interpolated conditions, IND: nominal order, DEP: nominal order
        {
            // Create a tabulated atmosphere object.
            std::string tabulatedAtmosphereFile = input_output::getAtmosphereTablesPath( ) +
                    "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
            TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFile );

            // Declare input and output
            altitudeInput = 10.05e3;
            densityInput = 0.4110;
            pressureInput = 26299.0;
            temperatureInput = 222.9350;

            densityOutput = tabulatedAtmosphere.getDensity( altitudeInput );
            pressureOutput = tabulatedAtmosphere.getPressure( altitudeInput );
            temperatureOutput = tabulatedAtmosphere.getTemperature( altitudeInput );
            break;
        }
        case 2: // 1-D: corner conditions, IND: shuffled order, DEP: shuffled order
        {
            // Create a tabulated atmosphere object.
            std::map< int, std::string > tabulatedAtmosphereFiles = { { 0, input_output::getAtmosphereTablesPath( ) +
                                                                        "MCDMeanAtmosphere.dat" } };
            std::vector< AtmosphereDependentVariables > dependentVariables = {
                pressure_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere, gas_constant_dependent_atmosphere,
                temperature_dependent_atmosphere, density_dependent_atmosphere };
            std::vector< AtmosphereIndependentVariables > independentVariables = { latitude_dependent_atmosphere };
            TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables );

            // Declare input and output
            latitudeInput = 5.0e4;
            densityInput = 1.4331262554e+00;
            pressureInput = 7.6178752157e-05;
            temperatureInput = 2.0951356051e+02;
            double gasConstantInput = 1.6104994141e+02;
            double ratioSpecificHeatInput = 2.3477457225e+00;
            double soundSpeedInput = 281.456889155598;

            densityOutput = tabulatedAtmosphere.getDensity( 0.0, 0.0, latitudeInput );
            pressureOutput = tabulatedAtmosphere.getPressure( 0.0, 0.0, latitudeInput );
            temperatureOutput = tabulatedAtmosphere.getTemperature( 0.0, 0.0, latitudeInput );
            double gasConstantOutput = tabulatedAtmosphere.getSpecificGasConstant( 0.0, 0.0, latitudeInput );
            double ratioSpecificHeatOutput = tabulatedAtmosphere.getRatioOfSpecificHeats( 0.0, 0.0, latitudeInput );
            double soundSpeedOutput = tabulatedAtmosphere.getSpeedOfSound( 0.0, 0.0, latitudeInput );

            std::cout << "Gas Constant Output: " << gasConstantOutput << " J/kg/K" << std::endl;
            std::cout << "Difference Gas Constant: " << gasConstantOutput - gasConstantInput << " J/kg/K" << std::endl << std::endl;
            std::cout << "Specific Heat Ratio Output: " << ratioSpecificHeatOutput << std::endl;
            std::cout << "Difference Specific Heat Ratio: " << ratioSpecificHeatOutput - ratioSpecificHeatInput << std::endl << std::endl;
            std::cout << "Speed Of Sound Output: " << soundSpeedOutput << " m/s" << std::endl;
            std::cout << "Difference Speed Of Sound: " << soundSpeedOutput - soundSpeedInput << " m/s" << std::endl << std::endl;
            break;
        }
        case 3: // 3-D: corner conditions, IND: shuffled order, DEP: nominal order
        {
            // Create a tabulated atmosphere object.
            std::map< int, std::string > tabulatedAtmosphereFiles;
            tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/density.dat";
            tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/pressure.dat";
            tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/temperature.dat";

            std::vector< AtmosphereDependentVariables > dependentVariables = {
                density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere };
            std::vector< AtmosphereIndependentVariables > independentVariables = {
                longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };

            TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables );

            // Declare input and output
            altitudeInput = 5.0e4;
            longitudeInput = convertDegreesToRadians( -180.0 );
            latitudeInput = convertDegreesToRadians( -90.0 );
            densityInput = 5.2805275e-05;
            pressureInput = 1.6627685;
            temperatureInput = 151.544;

            densityOutput = tabulatedAtmosphere.getDensity( altitudeInput, longitudeInput, latitudeInput );
            pressureOutput = tabulatedAtmosphere.getPressure( altitudeInput, longitudeInput, latitudeInput );
            temperatureOutput = tabulatedAtmosphere.getTemperature( altitudeInput, longitudeInput, latitudeInput );
            break;
        }
        case 4: // 3-D: non-interpolated conditions, IND: shuffled order, DEP: shuffled order
        {
            // Create a tabulated atmosphere object.
            std::map< int, std::string > tabulatedAtmosphereFiles;
            tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/temperature.dat";
            tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/density.dat";
            tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/pressure.dat";

            std::vector< AtmosphereDependentVariables > dependentVariables = {
                temperature_dependent_atmosphere, density_dependent_atmosphere, pressure_dependent_atmosphere };
            std::vector< AtmosphereIndependentVariables > independentVariables = {
                longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };

            TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables );

            // Declare input and output
            altitudeInput = 6.6202093e5;
            longitudeInput = convertDegreesToRadians( -1.685714e+02 );
            latitudeInput = convertDegreesToRadians( -7.851064e+01 );
            densityInput = 2.03357457566921e-15;
            pressureInput = 2.28335049368378e-09;
            temperatureInput = 204.242247282833;

            densityOutput = tabulatedAtmosphere.getDensity( altitudeInput, longitudeInput, latitudeInput );
            pressureOutput = tabulatedAtmosphere.getPressure( altitudeInput, longitudeInput, latitudeInput );
            temperatureOutput = tabulatedAtmosphere.getTemperature( altitudeInput, longitudeInput, latitudeInput );
            break;
        }
        case 5: // 3-D: interpolated conditions, IND: shuffled order, DEP: nominal order
        {
            // Create a tabulated atmosphere object.
            std::map< int, std::string > tabulatedAtmosphereFiles;
            tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/density.dat";
            tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/pressure.dat";
            tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/temperature.dat";
            tabulatedAtmosphereFiles[ 3 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/gasConstant.dat";
            tabulatedAtmosphereFiles[ 4 ] = input_output::getAtmosphereTablesPath( ) +
                    "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";

            std::vector< AtmosphereDependentVariables > dependentVariables = {
                density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere,
                gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere };
            std::vector< AtmosphereIndependentVariables > independentVariables = {
                longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };

            TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, independentVariables, dependentVariables );

            // Declare input and output
            altitudeInput = 236.9862e3;
            longitudeInput = convertDegreesToRadians( 72.98632 );
            latitudeInput = convertDegreesToRadians( -65.9762 );
            densityInput = 5.50931580592416e-13;
            pressureInput = 5.05339201226489e-08;
            temperatureInput = 174.82724294922;
            double gasConstantInput = 388.076687938436;
            double ratioSpecificHeatInput = 1.50561688237826;
            double soundSpeedInput = 319.610155078634;

            densityOutput = tabulatedAtmosphere.getDensity( altitudeInput, longitudeInput, latitudeInput );
            pressureOutput = tabulatedAtmosphere.getPressure( altitudeInput, longitudeInput, latitudeInput );
            temperatureOutput = tabulatedAtmosphere.getTemperature( altitudeInput, longitudeInput, latitudeInput );
            double gasConstantOutput = tabulatedAtmosphere.getSpecificGasConstant( altitudeInput, longitudeInput, latitudeInput );
            double ratioSpecificHeatOutput = tabulatedAtmosphere.getRatioOfSpecificHeats( altitudeInput, longitudeInput, latitudeInput );
            double soundSpeedOutput = tabulatedAtmosphere.getSpeedOfSound( altitudeInput, longitudeInput, latitudeInput );

            std::cout << "Gas Constant Output: " << gasConstantOutput << " J/kg/K" << std::endl;
            std::cout << "Difference Gas Constant: " << gasConstantOutput - gasConstantInput << " J/kg/K" << std::endl << std::endl;
            std::cout << "Specific Heat Ratio Output: " << ratioSpecificHeatOutput << std::endl;
            std::cout << "Difference Specific Heat Ratio: " << ratioSpecificHeatOutput - ratioSpecificHeatInput << std::endl << std::endl;
            std::cout << "Speed Of Sound Output: " << soundSpeedOutput << " m/s" << std::endl;
            std::cout << "Difference Speed Of Sound: " << soundSpeedOutput - soundSpeedInput << " m/s" << std::endl << std::endl;
            break;
        }
        }

        // Show tabulated output
        std::cout << "Density Output: " << densityOutput << " kg/m^3" << std::endl;
        std::cout << "Difference Density: " << densityOutput - densityInput << " kg/m^3" << std::endl << std::endl;
        std::cout << "Pressure Output: " << pressureOutput<< " Pa" << std::endl;
        std::cout << "Difference Pressure: " << pressureOutput- pressureInput << " Pa" << std::endl << std::endl;
        std::cout << "Temperature Output: " << temperatureOutput << " K" << std::endl;
        std::cout << "Difference Temperature: " << temperatureOutput - temperatureInput << " K" << std::endl;
        break;
    }
    case 1:
    {
        // Select file
        //      0: 4th dimension along rows
        //      1: 4th dimension along columns
        const int fileType = 0;
        std::string fileName;
        switch ( fileType )
        {
        case 0:
            fileName = tudat::input_output::getTudatRootPath( )
                    + "/Astrodynamics/Aerodynamics/UnitTests/dCDw4DTest1.txt";
            break;
        case 1:
            fileName = tudat::input_output::getTudatRootPath( )
                    + "/Astrodynamics/Aerodynamics/UnitTests/dCDw4DTest2.txt";
            break;
        }
        boost::multi_array< double, 4 > multiArrayFromFile =
                tudat::input_output::MultiArrayFileReader< 4 >::readMultiArray( fileName );

        // Extract data from file
        std::pair< boost::multi_array< double, 4 >, std::vector< std::vector< double > > > fileContents =
                tudat::input_output::MultiArrayFileReader< 4 >::readMultiArrayAndIndependentVariables( fileName );
        multiArrayFromFile = fileContents.first;
        std::vector< std::vector< double > > independentVariables = fileContents.second;

        // Test independent variables size
        std::cout << "# IV: " << independentVariables.size( ) << ", Exp: 4" << std::endl;
        std::cout << "1st IV: " << independentVariables.at( 0 ).size( ) << ", Exp: 11" << std::endl;
        std::cout << "2nd IV: " << independentVariables.at( 1 ).size( ) << ", Exp: 9" << std::endl;
        std::cout << "3rd IV: " << independentVariables.at( 2 ).size( ) << ", Exp: 5" << std::endl;
        std::cout << "4th IV: " << independentVariables.at( 3 ).size( ) << ", Exp: 2" << std::endl << std::endl;

        // Test multi-array size
        std::cout << "1st IV: " << multiArrayFromFile.shape( )[ 0 ] << ", Exp: 11" << std::endl;
        std::cout << "2nd IV: " << multiArrayFromFile.shape( )[ 1 ] << ", Exp: 9" << std::endl;
        std::cout << "3rd IV: " << multiArrayFromFile.shape( )[ 2 ] << ", Exp: 5" << std::endl;
        std::cout << "4th IV: " << multiArrayFromFile.shape( )[ 3 ] << ", Exp: 2" << std::endl << std::endl;

        // Test selected multi-array values in first 4th dimension
        std::cout << "1st Value: " << multiArrayFromFile[ 1 ][ 3 ][ 1 ][ 0 ] << ", Exp: 0" << std::endl;
        std::cout << "2nd Value: " << multiArrayFromFile[ 2 ][ 4 ][ 1 ][ 0 ] << ", Exp: -0.002" << std::endl;
        std::cout << "3rd Value: " << multiArrayFromFile[ 2 ][ 6 ][ 1 ][ 0 ] << ", Exp: -0.008" << std::endl;
        std::cout << "4th Value: " << multiArrayFromFile[ 1 ][ 3 ][ 3 ][ 0 ] << ", Exp: 0.0028" << std::endl << std::endl;

        // Test selected multi-array values in second 4th dimension
        std::cout << "1st Value: " << multiArrayFromFile[ 1 ][ 3 ][ 1 ][ 1 ] << ", Exp: 0" << std::endl;
        std::cout << "2nd Value: " << multiArrayFromFile[ 2 ][ 4 ][ 1 ][ 1 ] << ", Exp: -0.002" << std::endl;
        std::cout << "3rd Value: " << multiArrayFromFile[ 2 ][ 6 ][ 1 ][ 1 ] << ", Exp: -0.008" << std::endl;
        std::cout << "4th Value: " << multiArrayFromFile[ 1 ][ 3 ][ 3 ][ 1 ] << ", Exp: 0.0028" << std::endl;
        if ( fileType == 0 )
        {
            std::cout << "5th Value: " << multiArrayFromFile[ 9 ][ 6 ][ 4 ][ 1 ] << ", Exp: 0.2949" << std::endl;
        }
        else
        {
            std::cout << "5th Value: " << multiArrayFromFile[ 9 ][ 6 ][ 4 ][ 1 ] << ", Exp: 0.3949" << std::endl;
        }
        break;
    }
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

