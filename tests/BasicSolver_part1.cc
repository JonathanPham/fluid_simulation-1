#define _USE_MATH_DEFINES  // use math defines from std library

#include <iostream>
#include <math.h>
#include <cmath>

#include "Debug.hh"
#include "Array.hh"
#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "FluidSimulator.hh"

void setupTest1( Array & u, Array & v )
{
   u( 0, 0 ) = real( 1 );
   u( 1, 0 ) = real( 1 );
   u( 2, 0 ) = real( 1 );

   u( 0, 1 ) = real( 1 );
   u( 1, 1 ) = real( 2 );
   u( 2, 1 ) = real( 5 );

   u( 0, 2 ) = real( 1 );
   u( 1, 2 ) = real( 3 );
   u( 2, 2 ) = real( 9 );

   u( 0, 3 ) = real( 1 );
   u( 1, 3 ) = real( 5 );
   u( 2, 3 ) = real( 10 );

   v( 0, 0 ) = real( 1 );
   v( 1, 0 ) = real( 1 );
   v( 2, 0 ) = real( 1 );
   v( 3, 0 ) = real( 1 );

   v( 0, 1 ) = real( 1 );
   v( 1, 1 ) = real( 3 );
   v( 2, 1 ) = real( 6 );
   v( 3, 1 ) = real( 7 );

   v( 0, 2 ) = real( 1 );
   v( 1, 2 ) = real( 4 );
   v( 2, 2 ) = real( 8 );
   v( 3, 2 ) = real( 9 );
}

void checkTest1( Array & f, Array & g, int imax, int jmax )
{
   Array fSolution( imax + 1, jmax + 2 ), gSolution( imax + 2, jmax + 1 );

   fSolution.fill( real( 0 ) );

   // set values of F(i, j) from the analytical solution
   fSolution( 0, 1 ) = real( 1 );
   fSolution( 1, 1 ) = real( -0.745 );
   fSolution( 2, 1 ) = real( 5 );
   fSolution( 0, 2 ) = real( 1 );
   fSolution( 1, 2 ) = real( -2.195 );
   fSolution( 2, 2 ) = real( 9 );

   // compare computed solution of F with the analytical
   for( int j = 0; j <= ( jmax + 1 ); ++j )
   {
      for( int i = 0; i <= imax; ++i )
      {
         CHECK_REAL_EQUAL( f( i, j ), fSolution( i, j ) );
      }
   }

   gSolution.fill( real( 0 ) );
   // set values of F(i, j) from the analytical solution
   gSolution( 1, 0 ) = real( 1 );
   gSolution( 2, 0 ) = real( 1 );
   gSolution( 1, 1 ) = real( -0.05 );
   gSolution( 2, 1 ) = real( -8.76 );
   gSolution( 1, 2 ) = real( 4 );
   gSolution( 2, 2 ) = real( 8 );

   // compare computed solution of G with the analytical
   for( int j = 0; j <= jmax; ++j )
   {
      for( int i = 0; i <= ( imax + 1 ); ++i )
      {
         CHECK_REAL_EQUAL( g( i, j ), gSolution( i, j ) );
      }
   }
}

void setupTest2( Array & u, Array & v )
{
   u.fill( real( 1 ) );
   v.fill( real( -2 ) );
}

void checkTest2( Array & f, Array & g, Array & u, Array & v, int imax, int jmax, real dt, real gx, real gy )
{
   Array fSolution( imax + 1, jmax + 2 ), gSolution( imax + 2, jmax + 1 );

   fSolution.fill( real( 0 ) );

   // set up the solution array for f
   for( int j = 1; j <= jmax; ++j )  // set-up boundaries for f
   {
      fSolution( 0, j ) = u( 0, j );
      fSolution( imax, j ) = u( imax, j );
   }

   for( int j = 1; j <= jmax; ++j )
   {
      for( int i = 1; i < imax; ++i )
      {
         fSolution( i, j ) = u( i, j ) + dt * gx;
      }
   }

   // compare computed solution of F with the analytical
   for( int j = 0; j <= ( jmax + 1 ); ++j )
   {
      for( int i = 0; i <= imax; ++i )
      {
         CHECK_REAL_EQUAL( f( i, j ), fSolution( i, j ) );
      }
   }

   gSolution.fill( real( 0 ) );

   for( int i = 1; i <= imax; ++i )  // set-up boundaries for g
   {
      gSolution( i, 0 ) = v( i, 0 );
      gSolution( i, jmax ) = v( i, jmax );
   }

   for( int j = 1; j < jmax; ++j )
   {
      for( int i = 1; i <= imax; ++i )
      {
         gSolution( i, j ) = v( i, j ) + dt * gy;
      }
   }

   // compare computed solution of G with the analytical
   for( int j = 0; j <= jmax; ++j )
   {
      for( int i = 0; i <= ( imax + 1 ); ++i )
      {
         CHECK_REAL_EQUAL( g( i, j ), gSolution( i, j ) );
      }
   }
}

int main( int argc, char** argv )
{
   // check the number of parameters
   CHECK_MSG( argc >= 2, "You must supply a configuration file name!" );
   // config file name
   std::string confFile( argv[1] );
   // initialize conf. reader
   FileReader parameters;
   PROGRESS( "Reading input file '" << confFile << "'" );
   bool parametersRead = parameters.readFile( confFile );
   CHECK_MSG( parametersRead, "Failed to read input file '" + confFile + "'!" );
   // print cofiguration parameters
   parameters.printParameters();

   FluidSimulator simulator( parameters );

   Array & u = simulator.grid().u();
   Array & v = simulator.grid().v();
   Array & f = simulator.grid().f();
   Array & g = simulator.grid().g();

   // set up test case
   // assumptions for test case 1:
   //  gx = gy = 0
   //  gamma = 0.9
   //  Re = 10
   //  dt = 0.1
   //  imax = jmax = 2
   PROGRESS( "Setup test case 1" );
   f.fill( real( 0 ) );
   g.fill( real( 0 ) );
   setupTest1( u, v );
   // compute F and G
   PROGRESS( "Compute F and G for test case 1" );
   simulator.computeFG();
   // check solution
   PROGRESS( "Check solution for test case 1" );
   checkTest1( f, g, simulator.grid().imax(), simulator.grid().jmax() );
   PROGRESS( "Check for test case 1 succeeded" );

   // set up trivial test case
   // * note 1: for g_x = g_y = 0, the solution in F and G should be the corresponding values of u and v
   // * note 2: setting g_x and g_y to different values than 0 breaks test case 1
   PROGRESS( "Setup test case 2" );
   f.fill( real( 0 ) );
   g.fill( real( 0 ) );
   setupTest2( u, v );
   PROGRESS( "Compute F and G for test case 2" );
   simulator.computeFG();
   PROGRESS( "Check solution for test case 2" );
   checkTest2( f, g, u, v, simulator.grid().imax(), simulator.grid().jmax(), parameters.getRealParameter( "dt" ), parameters.getRealParameter( "gx" ),
               parameters.getRealParameter( "gy" ) );
   PROGRESS( "Check for test case 2 succeeded" );

   PROGRESS( "Program finished! Exiting..." );
   return 0;
}
