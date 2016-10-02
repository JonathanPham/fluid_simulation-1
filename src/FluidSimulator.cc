#include <limits>
#include <cmath>
#include "GrayScaleImage.hh"
#include "FluidSimulator.hh"

FluidSimulator::FluidSimulator( const FileReader & conf )
    : conf_( conf ), grid_( StaggeredGrid(conf) ), solver_( SORSolver(conf) ),
      vtk_( VTKWriter( grid_, conf.getStringParameter("name") ) )
{
    gamma_ = conf.getRealParameter( "gamma" );
    re_ = conf.getRealParameter("re");
    invRe_ = (real) 1 / re_;
    gx_ = conf.getRealParameter( "gx" );
    gy_ = conf.getRealParameter( "gy" );
    safetyFactor_ = conf.getRealParameter( "safetyfactor" );
    southBc_ = static_cast< BcType >( conf_.getIntParameter( "boundary_condition_s" ) );
    northBc_ = static_cast< BcType >( conf_.getIntParameter( "boundary_condition_n" ) );
    westBc_ = static_cast< BcType >( conf_.getIntParameter( "boundary_condition_w" ) );
    eastBc_ = static_cast< BcType >( conf_.getIntParameter( "boundary_condition_e" ) );

    checkParametersBC();
    createObstacle();
}

void FluidSimulator::simulate( real duration )
{
    real t = real(0),
        dt = conf_.getRealParameter( "dt" );

    unsigned long long n = 0;

    // get output interval and norm frequency parameters
    const int outIntv = conf_.getIntParameter("outputinterval"),
        normFreq = conf_.getIntParameter("normalizationfrequency");

    // assign initial values to u, v, p
    initFields();
    // write first VTK frame
    vtk_.write();
    while( t <= duration )
    {
        // simulate one timestep
        t += simulateOneTimestep( n, outIntv, normFreq, dt );
        if ( (n % outIntv) == 0 )
            PROGRESS("Simulation time: " << t);
        ++n;
    }
}

void FluidSimulator::simulateTimeStepCount( unsigned int nrOfTimeSteps )
{
    real t = real(0),
         dt = conf_.getRealParameter( "dt" );

    // get output interval and norm frequency parameters
    int outIntv = conf_.getIntParameter("outputinterval"),
        normFreq = conf_.getIntParameter("normalizationfrequency");
    // assign initial values to u, v, p
    initFields();
    // write first VTK frame
    vtk_.write();
    for( unsigned long long n = 0; n < nrOfTimeSteps; ++n )
    {
        t += simulateOneTimestep( n, outIntv, normFreq, dt );
        if ( (n % outIntv) == 0 )
            PROGRESS("Simulation time: " << t);
    }
    PROGRESS("Total simulation time = " << t);
}

inline real FluidSimulator::simulateOneTimestep( const int & n, const int & outIntv, const int & normFreq, const real & dt)
{
    // select dt
    real nextDt = determineNextDT( dt );
    // set boundary values for u and v
    refreshBoundaries();
    // compute F_n and G_n
    computeFG( nextDt );
    // compute RHS of the pressure equation
    composeRHS( nextDt );
    // solve the pressure equation
    solver_.solve( grid() );
    // compute u_(n+1) and v_(n+1)
    updateVelocities( nextDt );
    // write to VTK file
    if ( (n % outIntv) == 0)
        vtk_.write();
    // normalize pressure values around 0
    if ( (n % normFreq) == 0 )
        NormalizePressure();

    return nextDt;
}

/**
 * Compute F and G terms according to Eq. 3.36, 3.37
 */
void FluidSimulator::computeFG( real dt )
{
    adjustBoundariesFG();

    for( int j = 1; j <= grid_.jmax(); ++j )
    {
         for( int i = 1; i < grid_.imax(); ++i )
         {
             if (grid().isFluid(i, j))
                 grid().f(i, j) = grid().u(i, j) + dt * ( invRe_ * laplacianU(i, j) - duSqrDx(i, j) - duvDy(i, j) + gx_);
         }
    }

    for( int j = 1; j < grid_.jmax(); ++j )
    {
         for( int i = 1; i <= grid_.imax(); ++i )
         {
            if (grid().isFluid(i, j))
                grid().g(i, j) = grid().v(i, j) + dt * ( invRe_ * laplacianV(i, j) - dvSqrDy(i, j) - duvDx(i, j) + gy_ );
         }
    }
}

real FluidSimulator::laplacianU( const int & i, const int & j )
{
    const real dx = grid().dx(),
        dy = grid().dy();
    const real lx = (grid().u(i-1, j) - real(2) * grid().u(i, j) + grid().u(i+1, j)) / (dx * dx);
    const real ly = (grid().u(i, j-1) - real(2) * grid().u(i, j) + grid().u(i, j+1)) / (dy * dy);
    return lx + ly;
}

real FluidSimulator::duSqrDx( const int & i, const int & j)
{
    real res = real(0);
    const real iDx = real(1) / grid().dx();
    const real gDy = gamma_ / grid().dy();

    const real t0 = 0.5 * (grid().u(i, j) + grid().u(i+1, j));
    res += iDx * t0 * t0;
    const real t1 = 0.5 * (grid().u(i-1, j) + grid().u(i, j));
    res -= iDx * t1 * t1;
    res += gDy * fabs( t0 ) * 0.5 * (grid().u(i, j) - grid().u(i+1, j));
    res -= gDy * fabs( t1 ) * 0.5 * (grid().u(i-1, j) - grid().u(i, j));

    return res;
}

real FluidSimulator::duvDy( const int & i, const int & j)
{
    real res = real(0);
    const real iDy = real(1) / grid().dy(),
        gDy = gamma_ * iDy;

    const real t0 = 0.5 * (grid().v(i, j) + grid().v(i+1, j));
    const real t1 = 0.5 * (grid().v(i, j-1) + grid().v(i+1, j-1));

    res += iDy * t0 * 0.5 * (grid().u(i, j) + grid().u(i, j+1));
    res -= iDy * t1 * 0.5 * (grid().u(i, j-1) + grid().u(i, j));
    res += gDy * fabs(t0) * (grid().u(i, j) - grid().u(i, j+1));
    res -= gDy * fabs(t1) * (grid().u(i, j-1) - grid().u(i, j));

    return res;
}

real FluidSimulator::laplacianV( const int & i, const int & j )
{
    const real dx = grid().dx(),
        dy = grid().dy();
    const real lx = (grid().v(i-1, j) - real(2) * grid().v(i, j) + grid().v(i+1, j)) / (dx * dx);
    const real ly = (grid().v(i, j-1) - real(2) * grid().v(i, j) + grid().v(i, j+1)) / (dy * dy);
    return lx + ly;
}

real FluidSimulator::dvSqrDy( const int & i, const int & j)
{
    real res = real(0);

    const real iDy = real(1) / grid().dy(),
        gDy = gamma_ * iDy;

    const real t0 = 0.5 * (grid().v(i, j) + grid().v(i, j+1));
    res += iDy * t0 * t0;
    const real t1 = 0.5 * (grid().v(i, j-1) + grid().v(i, j));
    res -= iDy * t1 * t1;
    res += gDy * fabs(t0) * 0.5 * (grid().v(i, j) - grid().v(i, j+1));
    res -= gDy * fabs(t1) * 0.5 * (grid().v(i, j-1) - grid().v(i, j));

    return res;
}

real FluidSimulator::duvDx( const int & i, const int & j)
{
    real res = real(0);
    const real iDx = real(1) / grid().dx(),
        gDx = gamma_ * iDx;

    const real t0 = 0.5 * (grid().u(i, j) + grid().u(i, j+1));
    const real t1 = 0.5 * (grid().u(i-1, j) + grid().u(i-1, j+1));

    res += iDx * t0 * 0.5 * (grid().v(i, j) + grid().v(i+1, j));
    res -= iDx * t1 * 0.5 * (grid().v(i-1, j) + grid().v(i, j));
    res += gDx * fabs(t0) * 0.5 * (grid().v(i, j) - grid().v(i+1, j));
    res -= gDx * fabs(t1) * 0.5 * (grid().v(i-1, j) - grid().v(i, j));

    return res;
}

inline void FluidSimulator::adjustBoundariesFG()
{
    // set G(i, 0) = v(i, 0)
    copyBoundaryPointsFG( grid().g(), grid().v(), grid().imax(), grid().jmax(), SOUTH );
    // set G(i, jmax) = v(i, jmax)
    copyBoundaryPointsFG( grid().g(), grid().v(), grid().imax(), grid().jmax(), NORTH );
    // set F(0, j) = u(0, j)
    copyBoundaryPointsFG( grid().f(), grid().u(), grid().imax(), grid().jmax(), WEST );
    // set F(imax, j) = u(imax, j)
    copyBoundaryPointsFG( grid().f(), grid().u(), grid().imax(), grid().jmax(), EAST );
}

inline void FluidSimulator::copyBoundaryPointsFG( Array<real> & dst, Array<real> & src, int imax, int jmax, Direction dir )
{
    int iMin = 0, jMin = 0, iMax = imax, jMax = jmax;

    switch( dir )
    {
    case SOUTH:  // south boundary
        jMax = jMin;
        iMin = 1;
        break;
    case NORTH:  // north boundary
        jMin = jMax;
        iMin = 1;
        break;
    case WEST:  // west boundary
        iMax = iMin;
        jMin = 1;
        break;
    case EAST:  // east boundary
        iMin = iMax;
        jMin = 1;
        break;
    default:
        ERROR( "Boundary direction not supported! ");
    }

    for( int j = jMin; j <= jMax; ++j )
    {
        for( int i = iMin; i <= iMax; ++i )
        {
            // set boundary point
            dst( i, j ) = src( i, j );
        }
    }
}


void FluidSimulator::updateVelocities(real dt)
{
    const int imax = grid().imax();
    const int jmax = grid().jmax();

    // update velocity u
    const real dtDx = dt / grid_.dx();
    grid().umax( std::numeric_limits<real>::min() );
    for( int j = 1; j <= jmax; ++j )
    {
        for( int i = 1; i < imax; ++i ) {
            if ( !grid().isFluid(i, j) )
                continue;

            real u_new = grid().f( i, j ) - dtDx * ( grid().p( i+1, j ) - grid().p( i, j ) );

            real u_newAbs = fabs(u_new);
            if ( u_newAbs > grid_.umax() )
                grid_.umax( u_newAbs );

            grid().u( i, j ) = u_new;
        }
    }
    // update velocity v
    const real dtDy = dt / grid_.dy();
    grid_.vmax( real( std::numeric_limits<real>::min() ) );
    for( int j = 1; j < jmax; ++j )
    {
        for( int i = 1; i <= imax; ++i )
        {
            if ( !grid().isFluid(i, j) )
                continue;

            real v_new = grid().g( i, j ) - dtDy * ( grid().p( i, j+1 ) - grid().p( i, j ) );

            real v_newAbs = fabs( v_new );
            if ( v_newAbs > grid_.vmax() )
                grid_.vmax( v_newAbs );

            grid().v( i, j ) = v_new;
        }
    }
}

/**
 * Determine the next timestep dt according to equation (3.50) in order to satisfy
 * each of the three stability conditions given in (3.49). The stability conditions are:
 *
 * (i)   (2 * dt) / Re < 1 / ( 1 / (dx * dx) + 1 / (dy * dy) )
 * (ii)  |umax| * dt < dx
 * (iii) |vmax| * dt < dy
 */
real FluidSimulator::determineNextDT( real dt )
{
    real nextDT;

    // Use the specified dt if the safety factor tao is smaller
    // than zero. Then, one needs to be sure that for each timestep,
    // the stability conditions are satisfied.
    if ( safetyFactor_ < real(0) )
        return dt;

    const real dx = grid().dx(),
        dy = grid().dy(),
        uMax = grid().umax(),
        vMax = grid().vmax();

    // Evaluate stability condition (i)
    real t0 = (re_ / real(2)) * (1 / ((1 / (dx * dx)) + (1 / (dy * dy)) ));
    // Evaluate stability condition (ii)
    real t1 = dx / uMax;
    // Check min( (i), (ii) )
    nextDT = std::min<real>( t0, t1 );
    // Evaluate stability condition (iii)
    real t2 = dy / vMax;
    // Check min( min( (i), (ii) ), (iii) )
    nextDT = std::min<real>( nextDT, t2 );
    // return the new timestep dt
    return safetyFactor_ * nextDT;
}

void FluidSimulator::refreshBoundaries()
{
    real southVel = conf_.getRealParameter( "boundary_velocity_s" );
    refreshBoundary( southBc_, southVel, SOUTH );

    real northVel = conf_.getRealParameter( "boundary_velocity_n" );
    refreshBoundary( northBc_, northVel, NORTH );

    real westVel = conf_.getRealParameter( "boundary_velocity_w" );
    refreshBoundary( westBc_, westVel, WEST );

    real eastVel = conf_.getRealParameter( "boundary_velocity_e" );
    refreshBoundary( eastBc_, eastVel, EAST );
}

/**
 * Compose the right-hand side of the pressure equation, according to equation (3.38).
 */
void FluidSimulator::composeRHS(real dt)
{
    real invDt = real(1) / dt;
    real invDx = real(1) / grid_.dx();
    real invDy = real(1) / grid_.dy();

    for( auto j = 1; j <= grid().jmax(); ++j )
    {
        for( auto i = 1; i <= grid().imax(); ++i )
        {
            // equation (3.38)
            if (grid().isFluid(i, j))
                grid().rhs( i, j ) = invDt * ( ( grid().f(i, j) - grid().f(i-1, j) ) * invDx + ( grid().g(i, j) - grid().g(i, j-1) ) * invDy );
        }
    }
}

/**
 * Refresh the values for the boundary conditions on the domain walls
 */
inline void FluidSimulator::refreshBoundary( BcType boundaryCond, real dirVel, Direction dir)
{

    switch( boundaryCond )
    {
    case NOSLIP:
        // set no-slip boundary condition according to equations (3.21), (3.22)
        // and (3.23)
        setBcNoSlip( dirVel, dir );
        break;
    case SLIP:
        // set free-slip boundary condition according to equations (3.24) and (3.25)
        setBcFreeSlip( dirVel, dir );
        break;
    case OUTFLOW:
        if ( dirVel != real(0) )
            ERROR( "Velocity value cannot be set with outflow boundary condition!" );
        // set outflow boundary condition according to equation (3.26)
        setBcOutflow( dir );
        break;
    case INFLOW:
        // set inflow boundary condition similarly as in (3.22)
        setBcInflow( dirVel, dir );
        break;
    case PERIODIC:
        // set periodic boundary conditions
        setBcPeriodic( dir );
        break;
    default:
        ERROR( "Boundary condition not supported!" );
    }
}

/**
 * Set the no-slip boundary condition (Eq. 3.21, 3.22, 3.23)
 *
 * Set the no-slip boundary condition for a given direction +dir+ with velocity
 * +dirVel+ (the velocity is only relevant in case of a moving wall, e.g.
 * lid-driven cavity).
 *
 * The continuous velocities should vanish at the boundary to satisfy the no-slip
 * condition. For the values lying directly on the boundary, i.e. u at the east
 * and west direction, v at the south and north direction, we set the values to
 * zero.
 *
 * Since the vertical boundaries contain no v-values and the horizontal boundaries
 * contain no u-values, the zero boundary value is enforced in these cases by
 * averaging the values on either side of the boundary, e.g. on the west wall
 * for the v-values we set the point that would lie at the boundary v_r
 *   v_r = (v_a + v_i) / 2  <-- linear extrapolation
 * v_i is a fluid cell inside the domain and v_a is outside the domain.
 */
inline void FluidSimulator::setBcNoSlip( real dirVel, Direction dir )
{
    Array<real> & u = grid().u();
    Array<real> & v = grid().v();

    int imax = grid().imax(),
        jmax = grid().jmax();

    switch(dir)
    {
    case NORTH:
        // set values for velocity u
        for( auto i = 1; i <= imax; ++i)
            u( i, jmax + 1 ) = real(2) * dirVel - u( i, jmax );
        // set values for velocity v
        for( auto i = 1; i <= imax; ++i)
            v( i, jmax ) = real(0);
        break;
    case SOUTH:
        // set values for velocity u
        for( auto i = 1; i <= imax; ++i)
            u( i, 0 ) = real(2) * dirVel - u( i, 1 );
        // set values for velocity v
        for( auto i = 1; i <= imax; ++i)
            v( i, 0 ) = real(0);
        break;
    case WEST:
        // set values for velocity u
        for( auto j = 1; j <= jmax; ++j)
            u( 0, j ) = real(0);
        // set values for velocity v
        for( auto j = 1; j <= jmax; ++j)
            v( 0, j ) = real(2) * dirVel - v( 1, j );
        break;
    case EAST:
        // set values for velocity u
        for( auto j = 1; j <= jmax; ++j )
            u( imax, j ) = real(0);
        // set values for velocity v
        for( auto j = 1; j <= jmax; ++j )
            v( imax + 1, j ) = real(2) * dirVel - v( imax, j );
        break;
    default:
        ERROR( "Boundary direction not supported!" );
    }
}

/**
 * Set free slip boundary condition (Eq. 3.24, 3.25)
 *
 * In case of a free-slip boundary condition, the velocity component normal to
 * the boundary should vanish along with the normal derivative of the velocity
 * component tangent to the boundary. In the rectangular domain discretized
 * using the staggered grid, the values of velocities normal to the boundary
 * lie directly on the boundary; hence, just as for the no-slip condition, we
 * may set:
 *
 * u(0, j) = 0, u(imax, j) = 0, j = 1, ..., jmax,
 * v(i, 0) = 0, v(i, jmax) = 0, i = 1, ..., imax.
 *
 * For the tangential velocities, the discretized derivative of the velocity
 * with respect to the vector normal to the boundary needs to be calculated
 * and set equal to zero.
 */
inline void FluidSimulator::setBcFreeSlip( real dirVel, Direction dir )
{
    Array<real> & u = grid().u();
    Array<real> & v = grid().v();

    int imax = grid().imax(),
        jmax = grid().jmax();

    real dx = grid().dx(),
         dy = grid().dy();

    switch(dir)
    {
    case NORTH:
        //
        for( auto i = 1; i <= imax; ++i)
            u( i, jmax+1 ) = u( i, jmax ) + dy * dirVel;
        //
        for( auto i = 1; i <= imax; ++i)
            v( i, jmax ) = real(0);
        break;
    case SOUTH:
        //
        for( auto i = 1; i <= imax; ++i )
            u( i, 0 ) = u( i, 1 ) - dy * dirVel;

        //
        for( auto i = 1; i <= imax; ++i )
            v( i, 0 ) = real(0);
        break;
    case WEST:
        //
        for( auto j = 1; j <= jmax; ++j )
            u( 0, j ) = real(0);

        //
        for( auto j = 1; j <= jmax; ++j )
            v( 0, j ) = v( 1, j ) - dx * dirVel;
        break;
    case EAST:
        //
        for( auto j = 1; j <= jmax; ++j )
            u( imax, j ) = real(0);

        //
        for( auto j = 1; j <= jmax; ++j )
            v( imax+1, j ) = v( imax, j ) + dx * dirVel;
        break;
    default:
        ERROR( "Boundary direction not supported!" );
    }
}

/**
 * Set the outflow boundary condition (Eq. 3.26)
 *
 * In the outflow boundary condition the normal derivatives of both velocity
 * components are set to zero at the boundary, which means that the total velocity
 * does not change in the direction normal to the boundary. In the discrete case
 * this can be realized by setting velocity values at the boundary equal to their
 * neighbouring velocities inside the domain, i.e.,
 *
 *     u(0, j) = u(1, j),
 *     u(imax, j) = u(imax-1, j)
 *     v(0, j) = v(1, j),
 *     v(imax+1, j) = v(imax, j), for j = 1, ..., jmax
 *
 *     u(i, 0) = u(i, 1),
 *     u(i, jmax+1) = u(i, jmax)
 *     v(i, 0) = v(i, 1),
 *     v(i, jmax) = v(i, jmax-1), for i = 1, ..., imax
 */
inline void FluidSimulator::setBcOutflow( Direction dir )
{
    Array<real> & u = grid().u();
    Array<real> & v = grid().v();

    int imax = grid().imax(),
        jmax = grid().jmax();

    switch(dir)
    {
    case NORTH:
        // set values for velocity u
        for( auto i = 1; i <= imax; ++i )
            u( i , jmax + 1) = u( i, jmax );
        // set values for velocity v
        for( auto i = 1; i <= imax; ++i )
            v( i, jmax ) = v( i, jmax-1 );
        break;
    case SOUTH:
        // set values for velocity u
        for( auto i = 1; i <= imax; ++i )
            u( i, 0 ) = u( i, 1 );
        // set values for velocity v
        for( auto i = 1; i <= imax; ++i )
            v( i, 0 ) = v( i, 1 );
        break;
    case WEST:
        // set values for velocity u
        for( auto j = 1; j <= jmax; ++j )
            u( 0, j ) = u( 1, j );
        // set values for velocity v
        for( auto j = 1; j <= jmax; ++j )
            v( 0, j ) = v( 1, j );
        break;
    case EAST:
        // set values for velocity u
        for( auto j = 1; j <= jmax; ++j )
            u( imax, j ) = u( imax-1, j );
        // set values for velocity v
        for( auto j = 1; j <= jmax; ++j )
            v( imax+1, j ) = v( imax, j );
        break;
    default:
        ERROR( "Boundary direction not supported!" );
    }
}

/**
 * Set the inflow boundary condition
 *
 * On an inflow boundary the velocities are explicitly given; we impose this for
 * the velocities normal to the boundary (e.g., u on the left boundary) by directly
 * fixing the values on the boundary line. For the velocity components tangential
 * to the boundary (e.g., v at the left boundary), we achieve this by averaging
 * the values on either side of the boundary similarly as for the no-slip condition.
 **/
inline void FluidSimulator::setBcInflow( real dirVel, Direction dir )
{
    Array<real> & u = grid().u();
    Array<real> & v = grid().v();

    int imax = grid().imax(),
        jmax = grid().jmax();

    switch(dir)
    {
    case NORTH:
        // u is tangential to the upper boundary (north), we must
        // take the average of the inner point close to the boundary
        // u(i, jmax+1) and the point that lies outside the boundary
        // (similar to equation (3.22)).
        for( auto i = 1; i <= imax; ++i )
            u( i , jmax+1 ) = real(2) * dirVel - u( i, jmax );
        // v is normal to the upper boundary (north), we impose
        // values for u are equal the boundary velocity (dirVel)
        // in this direction.
        for( auto i = 1; i <= imax; ++i )
            v( i, jmax ) = dirVel;
        break;
    case SOUTH:
        // u is tangential to the upper boundary (north), we must
        // take the average of the inner point close to the boundary
        // u(i, 0) and the point that lies outside the boundary
        // (similar to equation (3.22)).
        for( auto i = 1; i <= imax; ++i )
            u( i, 0 ) = real(2) * dirVel - u( i, 1 );
        // v is normal to the upper boundary (north), we impose
        // values for u are equal the boundary velocity (dirVel)
        // in this direction.
        for( auto i = 1; i <= imax; ++i )
            v( i, 0 ) = dirVel;
        break;
    case WEST:
        // u is normal to the left boundary (west), we impose
        // values for u are equal the boundary velocity (dirVel)
        // in this direction.
        if ( conf_.getStringParameter("name") == "backstep" )
        {
            for( auto j = std::round(jmax / 2); j <= jmax; ++j )
                u( 0, j ) = dirVel;
        } else {
            for( auto j = 1; j <= jmax; ++j )
                u( 0, j ) = dirVel;
        }
        // v is tangential to the left boundary (west), we must
        // take the average of the inner point close to the boundary
        // v(1, j) and the point that lies outside the boundary (similar
        // to equation (3.22)).
        if ( !(conf_.getStringParameter("name") == "canal" || conf_.getStringParameter("name") == "backstep") )
        {
            for( auto j = 1; j <= jmax; ++j)
                v( 0, j ) = real(2) * dirVel - v( 1, j );
        }
        break;
    case EAST:
        // u is normal to the left boundary (west), we impose
        // values for u are equal the boundary velocity (dirVel)
        // in this direction.
        for( auto j = 1; j <= jmax; ++j )
            u( imax, j ) = dirVel;
        // v is tangential to the left boundary (west), we must
        // take the average of the inner point close to the boundary
        // v(imax+1, j) and the point that lies outside the boundary
        // (similar to equation (3.22)).
        for( auto j = 1; j <= jmax; ++j )
            v( imax + 1, j ) = real(2) * dirVel - v( imax, j );
        break;
    default:
        ERROR( "Boundary direction not supported!" );
    }
}

/**
 * Set periodic boundary condition (3.27)
 **/
inline void FluidSimulator::setBcPeriodic( Direction dir )
{
    Array<real> & u = grid().u();
    Array<real> & v = grid().v();

    int imax = grid().imax(),
        jmax = grid().jmax();

    Array<real> & p = grid().p();

    if ( dir == NORTH || dir == SOUTH )  // periodicity in x-direction
    {
        for( auto j = 1; j <= jmax; ++j )
        {
            u( 0, j ) = u( imax-1, j );
            u( 1, j ) = u( imax, j );
        }

        for( auto j = 1; j <= jmax; ++j )
        {
            v( 0, j ) = v(imax-1, j);
            v( 1, j ) = v(imax, j);
            v( 2, j ) = v(imax + 1, j);
        }

        for( auto j = 1; j <= jmax; ++j )
            p( 1, j ) = p( imax, j );
    }
    else if ( dir == EAST || dir == WEST )  // periodicity in y-direction
    {
        for( auto i = 1; i <= imax; ++i )
        {
            u( i, 0 ) = u(i, jmax-1);
            u( i, 1 ) = u(i, jmax);
            u( i, 2 ) = u(i, jmax+1);
        }

        for( auto i = 1; i <= imax; ++i )
        {
            v( i, 0 ) = v(i, jmax-1);
            v( i, 1 ) = v(i, jmax);
        }

        for( auto i = 1; i <= imax; ++i )
        {
            p( i, 1 ) = p(i, jmax);
        }
    }
    else  // should not be reached
    {
        ERROR( "Boundary direction not supported!" );
    }
}

void FluidSimulator::NormalizePressure()
{
    Array<real> & p = grid().p();
    int imax = grid_.imax(),
        jmax = grid_.jmax();

    const real normalization = real(1) / static_cast<real>(grid().getNumFluid());

    for( auto j = 1; j <= jmax; ++j )
    {
        for( auto i = 1; i <= imax; ++i )
        {
            p( i, j ) = p( i, j ) * normalization;
        }
    }
}

void FluidSimulator::initFields()
{
    grid().p().fill( conf_.getRealParameter( "p_init" ) );
    grid().u().fill( conf_.getRealParameter( "u_init" ) );
    grid().v().fill( conf_.getRealParameter( "v_init" ) );
    grid().resetFlagField();
}

void FluidSimulator::checkParametersBC()
{
    // check if both south and north boundaries are periodic
    if ( ((southBc_ == PERIODIC) || (northBc_ == PERIODIC)) && ( southBc_ != northBc_ ) )
    {
        ERROR("For periodicity in the y-direction, the north and south boundaries must be set to be periodic!");
    }
    // check if both west and east boundaries are periodic
    if ( ((westBc_ == PERIODIC) || (eastBc_ == PERIODIC)) && ( westBc_ != eastBc_ ) )
    {
        ERROR("For periodicity in the x-direction, the west and east boundaries must be set to be periodic!");
    }
}

void FluidSimulator::readDomainFromImg( const std::string & pngFilename )
{
    Array<bool> & flagField = grid().flagField();
    // read image from file
    GrayScaleImage img = GrayScaleImage( pngFilename );
    // resize image to fit domain
    img = img.getResizedImage( flagField.getSize(0), flagField.getSize(1) );

    for( int j = 0; j < flagField.getSize(1); ++j )
    {
        for( int i = 0; i < flagField.getSize(0); ++i )
            flagField( i, j ) = img(i, j) > real(0);
    }
}

void FluidSimulator::writeDomainToImg( const std::string & pngFilename )
{
    Array<bool> & flagField = grid().flagField();

    GrayScaleImage img = GrayScaleImage( flagField.getSize(0), flagField.getSize(1) );

    for( int j = 0; j < flagField.getSize(1); ++j )
    {
        for( int i = 0; i < flagField.getSize(0); ++i )
            img( i, j, static_cast<unsigned char>( (flagField(i, j)) ? 255 : 0 ) );
    }

    img.save( pngFilename );
}

void FluidSimulator::createObstacle()
{
    if ( conf_.getIntParameter("obstacletype") == int(RECTANGLE) ) {
        grid().createRectangle( conf_.getRealParameter("rectanglex1"), conf_.getRealParameter("rectangley1"),
          conf_.getRealParameter("rectanglex2"), conf_.getRealParameter("rectangley2") );
    } else if ( conf_.getIntParameter("obstacletype") == int(CIRCLE) ) {
        grid().createCircle( conf_.getRealParameter("circlex"), conf_.getRealParameter("circley"),
          conf_.getRealParameter("circler") );
    }
}
