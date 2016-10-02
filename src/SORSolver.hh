#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "Types.hh"
#include "StaggeredGrid.hh"


class SORSolver
{
public:
   // Constructor to manually create SORSolver
   SORSolver ( int itermax, real eps, real omega, int checkFrequency = int(1) );

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver ( const FileReader& configuration );

   // solve the pressure equation on the staggered grid
   bool solve( StaggeredGrid& grid );

private:
   int itermax_;
   int checkFrequency_;
   real eps_;
   real omega_;

   inline void SORIteration( StaggeredGrid & grid, real invDxSqr, real invDySqr, real invCenterWeight );
   inline real ResidualNorm( StaggeredGrid & grid, real invDxSqr, real invDySqr, real normalization);
   inline void setNeumannBC( StaggeredGrid & grid );
   inline void copyBoundaryPoints( StaggeredGrid & grid, Direction dir );
};






#endif //SOR_SOLVER_HH
