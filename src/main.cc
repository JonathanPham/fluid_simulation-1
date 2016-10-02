#define _USE_MATH_DEFINES  // use math defines from std library

#include <iostream>
#include <math.h>
#include <cmath>

#include "Debug.hh"
#include "Types.hh"
#include "Array.hh"
#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "FluidSimulator.hh"


int main( int argc, char** argv )
{
    // check the number of parameters
    CHECK_MSG( argc >= 2, "You must supply a configuration file name!" );
    // config file name
    std::string confFile( argv[1] );
    // initialize conf. reader
    FileReader parameters;
    //initParams( parameters );

    PROGRESS( "Reading input file '" << confFile << "'" );
    bool parametersRead = parameters.readFile( confFile );
    CHECK_MSG( parametersRead, "Failed to read input file '" + confFile + "'!" );

    parameters.printParameters();

    // set up fluid simulator
    FluidSimulator simulator( parameters );
    simulator.simulate( real(10) );

    return 0;
}
