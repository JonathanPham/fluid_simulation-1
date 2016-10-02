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

void initParams( FileReader & configuration )
{
    const std::map< std::string, real > realParameters = {
        /* Domain size */
        { "xlength", real(1) }, { "ylength", real(1) },
        /* Reynolds number */
        { "re", real(1) },
        /* Timestep size and safety factor (tao) */
        { "dt", real(0.01) }, { "safetyfactor", real(1) },
        /* External forces */
        { "gx", real(0) }, { "gy", real(0) },
        /* Velocities at the boundaries */
        { "boundary_velocity_n", real(0) }, { "boundary_velocity_s", real(0) },
        { "boundary_velocity_w", real(0) }, { "boundary_velocity_e", real(0) },
        /* Initial values for U, V, P */
        { "u_init", real(0) }, { "v_init", real(0) }, { "p_init", real(0) },
        /* Epsilon, Omega (SOR relaxation parameter), Gamma */
        { "eps", real(0) }, { "omg", real(0) }, { "gamma", real(0.1) }
    };

    const std::map< std::string, int > intParameters = {
        /* Domain mesh size in x, y directions */
        { "imax", int(1)  }, { "jmax", int(1) },
        /* Nr of timesteps */
        { "timesteps", int(100) },
        /* Boundary conditions for all directions */
        { "boundary_condition_n", NOSLIP }, { "boundary_condition_s", NOSLIP },
        { "boundary_condition_w", NOSLIP }, { "boundary_condition_e", NOSLIP },
        /* Number of iterations */
        { "itermax", int(1) },
        /* Residual check frequency */
        { "checkfrequency", int(1) },
        /* Pressure normalization frequency */
        { "normalizationfrequency", int(1) },
        /* VTK output frequency */
        { "outputinterval", int(1) }
    };

    const std::map< std::string, std::string > strParameters = {
        { "name", "" }
    };

    for( auto it = realParameters.begin(); it != realParameters.end(); ++it )
        configuration.registerRealParameter( it->first, it->second );

    for( auto it = intParameters.begin(); it != intParameters.end(); ++it )
        configuration.registerIntParameter( it->first, it->second );

    for( auto it = strParameters.begin(); it != strParameters.end(); ++it )
        configuration.registerStringParameter( it->first, it->second );
}

int main( int argc, char** argv )
{
    // check the number of parameters
    CHECK_MSG( argc >= 2, "You must supply a configuration file name!" );
    // config file name
    std::string confFile( argv[1] );
    // initialize conf. reader
    FileReader parameters;
    initParams( parameters );

    PROGRESS( "Reading input file '" << confFile << "'" );
    bool parametersRead = parameters.readFile( confFile );
    CHECK_MSG( parametersRead, "Failed to read input file '" + confFile + "'!" );

    parameters.printParameters();

    // set up fluid simulator
    FluidSimulator simulator( parameters );
    //simulator.simulateTimeStepCount( parameters.getIntParameter("timesteps") );
    simulator.simulate( real( 25 ) );
}
