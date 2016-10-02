#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH


#include "Types.hh"
#include "Array.hh"
#include "FileReader.hh"


//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
class StaggeredGrid
{
public:
    // Constructors to manually create staggered grid
    StaggeredGrid( int imax, int jmax, real dx, real dy, real umax = real(0), real vmax = real(0) );

    // Constructor to create a staggered grid from a parsed configuration file
    StaggeredGrid ( const FileReader & configuration );

    inline real & p( const int x, const int y )   { return p_(x, y); }
    inline real & rhs( const int x, const int y ) { return rhs_(x-1, y-1); }
    inline real & u( const int x, const int y ) { return u_(x, y); }
    inline real & v( const int x, const int y ) { return v_(x, y); }
    inline real & f( const int x, const int y ) { return f_(x, y); }
    inline real & g( const int x, const int y ) { return g_(x, y); }
    inline bool & flagField(const int x, const int y ) { return ff_(x, y); }
    real u( const int x, const int y, Direction dir );
    real v( const int x, const int y, Direction dir );
    real p( const int x, const int y, Direction dir );
    // Getters / Setters for member variables
    inline Array<real> & p()    { return p_;    }
    inline Array<real> & rhs()  { return rhs_;  }
    inline Array<real> & u()    { return u_; }
    inline Array<real> & v()    { return v_; }
    inline Array<real> & f()    { return f_; }
    inline Array<real> & g()    { return g_; }
    inline Array<bool> & flagField()   { return ff_; }

    inline const Array<real> & p()   const { return p_;   }
    inline const Array<real> & rhs() const { return rhs_; }
    inline const Array<real> & u()   const { return u_;   }
    inline const Array<real> & v()   const { return v_;   }
    inline const Array<real> & f()   const { return f_;   }
    inline const Array<real> & g()   const { return g_;   }
    inline const Array<bool> & ff()  const { return ff_; }

    inline int imax() const { return imax_; }
    inline int jmax() const { return jmax_; }

    inline real dx() const { return dx_; }
    inline real dy() const { return dy_; }

    inline real umax() const { return umax_; }
    inline real vmax() const { return vmax_; }
    inline void umax( real umax ) { CHECK( umax != real(0) ); umax_ = umax; }
    inline void vmax( real vmax ) { CHECK( vmax != real(0) ); vmax_ = vmax; }

    //< obstacle handling
    inline bool isFluid( const int x, const int y );
    inline int getNumFluid() { return numFluid_; }
    inline void setCellToObstacle( const int x, const int y );
    inline void resetFlagField() { ff_.fill(true); numFluid_ = imax_ * jmax_; }

    //< obstacle creation for basic geometries
    void createRectangle( real x1, real y1, real x2, real y2 );
    void createCircle( real x, real y, real r );

protected:
    Array<real> p_;   //< pressure field
    Array<real> rhs_; //< right hand side of the pressure equation
    Array<real> u_;   //< velocity field u
    Array<real> v_;   //< velocity field v
    Array<real> f_;
    Array<real> g_;
    Array<bool> ff_;

    int imax_;  //< imax
    int jmax_;  //< jmax
    int numFluid_; //< number of fluid cells

    real dx_;   //< distance between two grid points in x direction
    real dy_;   //< distance between two grid points in y direction

    real umax_;  //< maximum absolute value of u velocity
    real vmax_;  //< maximum absolute value of v velocity

    real xlength_; //< domain length in x direction
    real ylength_; //< domain length in y direction

    inline void createArrays(int xMax, int yMax);
};

bool StaggeredGrid::isFluid( const int x, const int y )
{
     ASSERT((x < ff_.getSize(0)) && (y < ff_.getSize(1)));
     return ff_(x, y);
}

void StaggeredGrid::setCellToObstacle( const int x, const int y )
{
    ASSERT((x < ff_.getSize(0)) && (y < ff_.getSize(1)));
    ff_(x, y) = false;
    --numFluid_;
}

#endif //STAGGERED_GRID_HH
