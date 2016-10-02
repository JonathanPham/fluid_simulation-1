#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "Types.hh"
#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "VTKWriter.hh"

class FluidSimulator
{
public:
    FluidSimulator( const FileReader & conf );

    /// Simulates a given time-length
    void simulate             ( real duration              );
    void simulateTimeStepCount( unsigned int nrOfTimeSteps );

    // Getter functions for the internally stored StaggeredGrid
    inline StaggeredGrid & grid() { return grid_; }
    inline const StaggeredGrid & grid() const { return grid_; }

    // read domain obstacle from image file
    void readDomainFromImg( const std::string & pngFilename );
    // write domain obstacle to image file
    void writeDomainToImg( const std::string & pngFilename );

private:
    const FileReader & conf_;
    StaggeredGrid grid_;
    SORSolver solver_;
    VTKWriter vtk_;
    real gamma_;
    real re_;
    real invRe_;
    real gx_;
    real gy_;
    real safetyFactor_;
    BcType southBc_;
    BcType northBc_;
    BcType westBc_;
    BcType eastBc_;

    real simulateOneTimestep( const int & n, const int & outIntv, const int & normFreq, const real & dt);

    void computeFG( real dt );

    real laplacianU( const int & i, const int & j );

    real duSqrDx( const int & i, const int & j);

    real duvDy( const int & i, const int & j);

    real laplacianV( const int & i, const int & j );

    real dvSqrDy( const int & i, const int & j);

    real duvDx( const int & i, const int & j);

    void adjustBoundariesFG();

    void copyBoundaryPointsFG( Array<real> & dst, Array<real> & src, int imax, int jmax, Direction dir );

    void updateVelocities(real dt);

    real determineNextDT( real dt );

    void refreshBoundaries();

    void composeRHS(real dt);

    // refresh boundary for a given direction dir
    inline void refreshBoundary( BcType boundaryCond, real dirVel, Direction dir);
    // set no-slip BC for a given direction dir
    inline void setBcNoSlip( real dirVel, Direction dir );
    // set free-slip BC for a given direction dir
    inline void setBcFreeSlip( real dirVel, Direction dir );
    // set outflow BC for a given direction dir
    inline void setBcOutflow( Direction dir );
    // set inflow BC for a given direction dir
    inline void setBcInflow( real dirVel, Direction dir );
    // set periodic BC for a given direction dir
    inline void setBcPeriodic( Direction dir );
    // normalize pressure values
    void NormalizePressure();
    // initialize velocity fields u, v and pressure p
    void initFields();
    // check boundary condition parameters
    void checkParametersBC();
    // create obstacle inside the domain
    void createObstacle();
};



#endif
