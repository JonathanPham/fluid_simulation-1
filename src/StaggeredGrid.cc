#include <cmath>
#include <functional>
#include "Debug.hh"
#include "Array.hh"
#include "StaggeredGrid.hh"

StaggeredGrid::StaggeredGrid( int imax, int jmax, real dx, real dy, real umax, real vmax )
    : imax_(imax), jmax_(jmax), dx_(dx), dy_(dy), umax_(umax), vmax_(vmax),
      xlength_(dx * imax), ylength_(dy * jmax)
{
    ASSERT( imax > 0 );
    ASSERT( jmax > 0 );
    createArrays( imax_, jmax_ );
}

StaggeredGrid::StaggeredGrid ( const FileReader & configuration )
{
   // calculate d{x,y}
   xlength_ = configuration.getRealParameter("xlength");
   ylength_ = configuration.getRealParameter("ylength");

   imax_ = configuration.getIntParameter("imax");
   jmax_ = configuration.getIntParameter("jmax");

   ASSERT( imax_ > 0 );
   ASSERT( jmax_ > 0 );

   dx_ = xlength_ / static_cast< real >(imax_);
   dy_ = ylength_ / static_cast< real >(jmax_);

   umax_ = std::fabs( configuration.getRealParameter( "u_init" ) );
   vmax_ = std::fabs( configuration.getRealParameter( "v_init" ) );

   createArrays( imax_, jmax_ );
}

inline void StaggeredGrid::createArrays(int xMax, int yMax)
{
   rhs_ = Array< real >( xMax    , yMax     );
   p_   = Array< real >( xMax + 2, yMax + 2 );
   u_   = Array< real >( xMax + 1, yMax + 2 );
   v_   = Array< real >( xMax + 2, yMax + 1 );
   f_   = Array< real >( xMax + 1, yMax + 2 );
   g_   = Array< real >( xMax + 2, yMax + 1 );
   ff_  = Array< bool >( xMax + 2, yMax + 2 );
   resetFlagField();
}

/**
 * Wrapped access with direction to velocity field u
 */
real StaggeredGrid::u( const int x, const int y, Direction dir )
{
    if (isFluid(x, y) && isFluid(x+1, y))
        return u(x, y);

    switch(dir)
    {
    case NORTH:
        if (!isFluid(x, y) && !isFluid(x+1, y)) {
           // both SOLID, i.e., u(x, y) inside obstacle
            return -u(x, y-1);
        } else {
            // velocity u(x, y) is in a boundary, therefore we set it to 0 because of the no-slip cond.
            return real(0);
        }
        break;
    case SOUTH:
        if (!isFluid(x, y) && !isFluid(x+1, y)) {
            // both SOLID, i.e., u(x, y) inside obstacle
            return -u(x, y+1);
        } else {
            // velocity u(x, y) is in a boundary, therefore we set it to 0 because of the no-slip cond.
             return real(0);
        }
        break;
    case WEST:
        // cell u(x+1, y) cannot be solid, therefore it can either be that
        // cell (x, y) is an obstacle, then one must return 0. Otherwise,
        // the value is provided in the trivial case (both cells are fluid
        // cells).
        return real(0);
        break;
    case EAST:
        // if cell (x+1, y) is solid, then there will be no exchange of mass/energy
        // through the obstacle cell, since we're using a no-slip condition for the
        // obstacle.
        return real(0);
        break;
    default:  // should not be reached
        ERROR("Invalid stencil direction!");
        return real(0); //< prevent compiler warning (shouldn't be reached)

    }
}

/**
 * Wrapped access to velocity field v
 */
real StaggeredGrid::v( const int x, const int y, Direction dir )
{
    // trivial case, both velocities v(x, y) and v(x, y+1) are fluid cells.
    if ( isFluid(x, y) && isFluid(x, y+1) )
        return v( x, y );

    switch( dir )
    {
    case NORTH:
        // velocity v(x, y) can be either at a boundary (i.e. cell (x, y) is
        // fluid and (x, y+1) is solid) or (x, y) is solid, so there is no
        // exchange of fluid through the obstacle cell
        return real(0);
        break;
    case SOUTH:
        // velocity v(x, y+1) cannot be an obstacle cell, so v(x, y) can either
        // be a fluid cell, and this case is treated in the trivial case, or
        // an obstacle cell, in which case the velocity is at a boundary so
        // we must set it to 0.
        return real(0);
        break;
    case WEST:
        //
        if (!isFluid(x, y) && !isFluid(x, y+1)) { // velocity v(x, y) is inside the obstacle
            return -v(x+1, y);
        } else {
            return real(0);
        }
        break;
    case EAST:
        if (!isFluid(x, y) && !isFluid(x, y+1)) {  // velocity v(x, y) is inside the obstacle
            return -v(x-1, y);
        } else {  // velocity v(x, y) is at a boundary, therefore we set it to 0
            return real(0);
        }
        break;
    default:  // should not be reached
        ERROR("Invalid stencil direction!");
        return real(0); //< prevent compiler warning (shouldn't be reached)
    }
}

/**
 * Wrapped access to pressure array
 */
real StaggeredGrid::p( const int x, const int y, Direction dir )
{
    ASSERT_MSG( x < p_.getSize(0), "Access in x-direction is out-of-bounds! Requested index = " << x << " | Length in x-direction = " << p_.getSize(0) );
    ASSERT_MSG( y < p_.getSize(1), "Access in y-direction is out-of-bounds! Requested index = " << y << " | Length in y-direction = " << p_.getSize(1) );
    // check if cell (x, y) is fluid cell
    if ( isFluid( x, y ) )
        return p( x, y );
    // if cell(x, y) is not a fluid cell, compute the offsets for the value
    // of the pressure at the obstacle, which is *almost* equivalent to copying
    // the pressure values inside the obstacle cells
    switch( dir )
    {
    case NORTH: return p(x, y - 1); break;
    case EAST:  return p(x - 1, y); break;
    case SOUTH: return p(x, y + 1); break;
    case WEST:  return p(x + 1, y); break;
    default: ERROR("Direction not supported!"); return real(0);
    }
}

/**
 * Create rectangular obstacle inside the domain
 *
 * The rectangle is defined by two points (x1, y1) and (x2, y2)
 *
 * (x1, y2) ->  +---------------------+ <- (x2, y2)
 *              |                     |
 *              |                     |
 * (x1, y1) ->  +---------------------+ <- (x2, y1)
 **/
void StaggeredGrid::createRectangle( real x1, real y1, real x2, real y2 )
{
    CHECK_MSG( (x1 >= real(0)) && (x1 <= xlength_), "RectangleX1 is outside the domain!" );
    CHECK_MSG( (y1 >= real(0)) && (y1 <= ylength_), "RectangleY1 is outside the domain!" );
    CHECK_MSG( (x2 >= real(0)) && (x2 <= xlength_), "RectangleX2 is outside the domain!" );
    CHECK_MSG( (y2 >= real(0)) && (y2 <= ylength_), "RectangleY2 is outside the domain!" );
    CHECK_MSG( ( y1 < y2 ) && (x1 < x2), "Points (" << x1 << ", " << y1 << ") and ("
        << x2 << ", " << y2 << ") does not form a rectangle" );

    PROGRESS("Create rectangular obstacle with coordinates:\n\t(x1, y1) = (" << x1 << ", " << y1 << ")\n\t(x2, y2) = (" << x2 << ", " << y2 << ")");

    const int xEnd = std::floor( x2 / dx_ );
    const int yEnd = std::floor( y2 / dy_ );

    for( int j = std::floor( y1 / dy_ ); j <= yEnd; ++j )
    {
        for( int i = std::floor( x1 / dx_ ); i <= xEnd; ++i )
        {
            setCellToObstacle( i, j );
        }
    }
}


/**
 * Create a circle obstacle inside the domain
 *
 * The function assumes that coordinates x and y specify the center point of the
 * circle.
 */
void StaggeredGrid::createCircle( real x, real y, real r )
{
    CHECK_MSG( ( x >= 0 ) && ( x <= xlength_ ) && (y >= 0) && (y <= ylength_), "Point (" << x << ", " << y << ") is outside the domain!" );
    CHECK_MSG( ((x-r) >= 0) && ((x-r) <= xlength_) && ((x+r) >= 0) && ((x+r) <= xlength_), "Circle is outside the domain!" );
    CHECK_MSG( ((y-r) >= 0) && ((y-r) <= ylength_) && ((y+r) >= 0) && ((y+r) <= ylength_), "Circle is outside the domain!" );

    PROGRESS("Create circle obstacle with coordinates:\n\t(x, y) = (" << x << ", " << y << ")\n\tRadius: " << r);

    const int yEnd = std::floor( ( y + r ) / dy_ );
    const int xEnd = std::floor( ( x + r ) / dx_ );

    // lambda function to test if point is inside the circle
    std::function< bool( real xPoint, real yPoint ) > evalPoint = [x, y, r]( real xPoint, real yPoint )
        { return (((xPoint-x)*(xPoint-x) + (yPoint-y)*(yPoint-y)) <= (r * r) ); };

    // assumption: point (x, y) in the middle of the cell
    for( int j = std::floor( ( y - r ) / dy_ ); j <= yEnd; ++j )
    {
        for( int i = std::floor( ( x - r ) / dx_ ); i <= xEnd; ++i )
        {
            // test if all for points of cell (x, y) are inside the circle,
            // if yes then set cell to obstacle
            if ( evalPoint( i - 0.5 * dx_, j - 0.5 * dy_ ) && evalPoint( i + 0.5 * dx_, j - 0.5 * dy_ ) &&
                 evalPoint( i - 0.5 * dx_, j + 0.5 * dy_ ) && evalPoint( i + 0.5 * dx_, j + 0.5 * dy_ ) )
            {
                setCellToObstacle( i, j );
            }
        }
    }
}
