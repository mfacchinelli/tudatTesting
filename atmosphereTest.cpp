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

    // Declare input and output
    double altitudeInput, longitudeInput, latitudeInput, densityInput, pressureInput, temperatureInput;
    double densityOutput, pressureOutput, temperatureOutput;

    // Select case
    //      0: 1-D: corner conditions, IND: nominal order, DEP: nominal order
    //      1: 1-D: interpolated conditions, IND: nominal order, DEP: nominal order
    //      2: 1-D: corner conditions, IND: nominal order, DEP: shuffled order
    //      3: 3-D: corner conditions, IND: shuffled order, DEP: nominal order
    //      4: 3-D: non-interpolated conditions, IND: shuffled order, DEP: shuffled order
    //      5: 3-D: interpolated conditions, IND: shuffled order, DEP: nominal order
    const int testCase = 5;
    std::cout << "Test Case < " << testCase << " > Selected" << std::endl;
    std::cout << "Tolerance is set at: " << tolerance << std::endl << std::endl;
    switch ( testCase )
    {
    case 0: // 1-D: 0 km altitude, IND: nominal order, DEP: nominal order
    {
        // Create a tabulated atmosphere object.
        std::map< int, std::string > tabulatedAtmosphereFiles = { { 0, input_output::getAtmosphereTablesPath( ) +
                                                                    "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" } };
        TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles );

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
        std::map< int, std::string > tabulatedAtmosphereFiles = { { 0, input_output::getAtmosphereTablesPath( ) +
                                                                    "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" } };
        TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles );

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
    case 2: // 1-D: 0 km altitude, IND: nominal order, DEP: shuffled order
    {
        // Create a tabulated atmosphere object.
        std::map< int, std::string > tabulatedAtmosphereFiles = { { 0, input_output::getAtmosphereTablesPath( ) +
                                                                    "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" } };
        std::vector< AtmosphereDependentVariables > dependentVariables = {
            temperature_dependent_atmosphere, density_dependent_atmosphere, pressure_dependent_atmosphere };
        TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, dependentVariables );

        // Declare input and output
        altitudeInput = 0.0;
        densityInput = 101320.0;
        pressureInput = 288.15;
        temperatureInput = 1.225;

        densityOutput = tabulatedAtmosphere.getDensity( altitudeInput );
        pressureOutput = tabulatedAtmosphere.getPressure( altitudeInput );
        temperatureOutput = tabulatedAtmosphere.getTemperature( altitudeInput );
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

        TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, dependentVariables, independentVariables );

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
                "MCDMeanAtmosphereTimeAverage/density.dat";
        tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) +
                "MCDMeanAtmosphereTimeAverage/pressure.dat";
        tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) +
                "MCDMeanAtmosphereTimeAverage/temperature.dat";

        std::vector< AtmosphereDependentVariables > dependentVariables = {
            temperature_dependent_atmosphere, density_dependent_atmosphere, pressure_dependent_atmosphere };
        std::vector< AtmosphereIndependentVariables > independentVariables = {
            longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };

        TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, dependentVariables, independentVariables );

        // Declare input and output
        altitudeInput = 3.739610e8;
        longitudeInput = convertDegreesToRadians( -1.685714e+02 );
        latitudeInput = convertDegreesToRadians( -7.851064e+01 );
        densityInput = 2.6315275e-12;
        pressureInput = 204.24225;
        temperatureInput = 1.486543e-18;

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

        std::vector< AtmosphereDependentVariables > dependentVariables = {
            density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere };
        std::vector< AtmosphereIndependentVariables > independentVariables = {
            longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };

        TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, dependentVariables, independentVariables );

        // Declare input and output
        altitudeInput = 236.9862e3;
        longitudeInput = convertDegreesToRadians( 72.98632 );
        latitudeInput = convertDegreesToRadians( -65.9762 );
        densityInput = 5.57022396263159e-13;
        pressureInput = 5.08715805228075e-08;
        temperatureInput = 174.826970929597;

        densityOutput = tabulatedAtmosphere.getDensity( altitudeInput, longitudeInput, latitudeInput );
        pressureOutput = tabulatedAtmosphere.getPressure( altitudeInput, longitudeInput, latitudeInput );
        temperatureOutput = tabulatedAtmosphere.getTemperature( altitudeInput, longitudeInput, latitudeInput );
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

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

