#include <cmath>
#include <omp.h>
#include "Types.hh"
#include "SORSolver.hh"

SORSolver::SORSolver( int itermax, real eps, real omega, int checkFrequency )
         : itermax_( itermax ), checkFrequency_(checkFrequency), eps_( eps ), omega_( omega )
{
}

SORSolver::SORSolver( const FileReader& configuration )
{
   itermax_ = configuration.getIntParameter( "itermax" );
   eps_ = configuration.getRealParameter( "eps" );
   omega_ = configuration.getRealParameter( "omg" );
   checkFrequency_ = configuration.getIntParameter("checkfrequency");
}

/**
 * Solve Poisson equation of the pressure according to Eq. (3.38)
 */
bool SORSolver::solve( StaggeredGrid & grid )
{
   // L-squared norm
   real l2Norm;
   // get grid parameters
   real dx = grid.dx();
   real dy = grid.dy();

   // pre-calculate constants to avoid divisions inside loops
   const real invDxSqr = 1.0 / (dx * dx);
   const real invDySqr = 1.0 / (dy * dy);
   const real invCenterWeight = 1.0 / ( (2.0 * invDxSqr) + (2.0 * invDySqr) );
   const real normalization = real(1) / static_cast<real>(grid.getNumFluid());

   for( int iter = 0; iter < itermax_; ++iter )
   {
      // (1) set Neumann boundary conditions for the pressure
      setNeumannBC( grid );
      // (2) perform one iteration of the SOR-method
      SORIteration( grid, invDxSqr, invDySqr, invCenterWeight );
      // (3) check residuum calculation frequency
      if ( (iter % checkFrequency_) == 0 )
      {
          // (4) set Neumann boundary conditions for the pressure
          // before calculating the residual
          setNeumannBC( grid );
          // (5) calculate residual norm
          l2Norm = ResidualNorm( grid, invDxSqr, invDySqr, normalization );
          // (6) check residual norm
          if( l2Norm < eps_ )
             return true;
      }
   }
   return false;
}

/**
 * Succesive overrelaxation (SOR) method according to Eq. (3.44)
 */
inline void SORSolver::SORIteration( StaggeredGrid & grid, real invDxSqr, real invDySqr, real invCenterWeight )
{
   const real om_omg = 1.0 - omega_;
   const real ct_omg = omega_ * invCenterWeight;
   const real cx_omg = ct_omg * invDxSqr;
   const real cy_omg = ct_omg * invDySqr;

   for( int rb = 0; rb < 2; ++rb)  // 0 = red sweep, 1 = black sweep
   {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
      for( int j = 1; j <= grid.jmax(); ++j)
      {
         // check if column is odd, if so check if red or black sweep
         int iBegin = ( (j & 1) == 0 ) ? ((rb == 0) ? 1 : 2) : ( (rb == 1) ? 2 : 1 );

         for( int i = iBegin; i <= grid.imax(); i += 2 )
         {
            grid.p(i, j) = om_omg * grid.p(i, j);
            grid.p(i, j) += cx_omg * grid.p(i-1, j, WEST);
            grid.p(i, j) += cx_omg * grid.p(i+1, j, EAST);
            grid.p(i, j) += cy_omg * grid.p(i, j-1, SOUTH);
            grid.p(i, j) += cy_omg * grid.p(i, j+1, NORTH);
            grid.p(i, j) += ct_omg * -grid.rhs(i, j);
         }
      }
   }
}

/**
 * Compute L-2 norm of the residual according to equation (3.46)
 */
inline real SORSolver::ResidualNorm( StaggeredGrid & grid, real invDxSqr, real invDySqr, real normalization)
{
   real sqrNorm = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) reduction(+:sqrNorm)
#endif
   for( int j = 1; j <= grid.jmax(); ++j )
   {
      for( int i = 1; i <= grid.imax(); ++i )
      {
         real res = grid.rhs(i, j);
         // compute \dfrac{\delta^{2} p}{\delta x^{2}}
         real t0 = invDxSqr * -grid.p(i-1, j, WEST);
         t0 += invDxSqr * 2.0 * grid.p(i, j);
         t0 += invDxSqr * -grid.p(i+1, j, EAST);
         res += t0;
         // compute \dfrac{\delta^{2} p}{\delta y^{2}}
         real t1 = invDySqr * -grid.p(i, j-1, SOUTH);
         t1 += invDySqr * 2.0 * grid.p(i, j);
         t1 += invDySqr * -grid.p(i, j+1, NORTH);
         res += t1;

         sqrNorm += res * res;
      }
   }
   // calculate normalized L-squared norm
   return std::sqrt( normalization * sqrNorm );
}

inline void SORSolver::setNeumannBC( StaggeredGrid & grid )
{
   copyBoundaryPoints( grid, SOUTH ); // south boundary
   copyBoundaryPoints( grid, NORTH ); // north boundary
   copyBoundaryPoints( grid, WEST ); // west boundary
   copyBoundaryPoints( grid, EAST ); // east boundary
}

inline void SORSolver::copyBoundaryPoints( StaggeredGrid & grid, Direction dir )
{
   int imin = 1,
       jmin = 1,
       imax = grid.imax(),
       jmax = grid.jmax(),
       ioff = 0,
       joff = 0;

   switch( dir )
   {
   case SOUTH:  // south boundary
      jmax = 1;
      joff = -1;
      break;
   case NORTH:  // north boundary
      jmin = jmax;
      joff = 1;
      break;
   case WEST:  // west boundary
      imax = imin;
      ioff = -1;
      break;
   case EAST:  // east boundary
      imin = imax;
      ioff = 1;
      break;
   default:
      ERROR( "Boundary direction not supported!" );
   }

   for( int j = jmin; j <= jmax; ++j )
   {
      for( int i = imin; i <= imax; ++i )
      {
         // set boundary point
         grid.p(i+ioff, j+joff) = grid.p(i, j, dir);
      }
   }
}
