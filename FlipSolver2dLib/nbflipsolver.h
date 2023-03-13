#ifndef NBFLIPSOLVER_H
#define NBFLIPSOLVER_H

#include "flipsolver2d.h"
#include "sdfgrid.h"
#include "staggeredvelocitygrid.h"

class NBFlipSolver : public FlipSolver
{
public:
    NBFlipSolver(int sizeI, int sizeJ, int extrapRadius = 1, bool vonNeumannNeighbors = false);

protected:
    void step() override;

    void advect() override;

    void afterTransfer() override;

    void particleVelocityToGrid() override;

    void reseedParticles() override;

    void firstFrameInit() override;

    void initialFluidSeed();

    void fluidSdfFromInitialFluid();

    void combineAdvectedGrids();

    void combineLevelset();

    void combineVelocityGrid();

    void extrapolateLevelsetInside(SdfGrid& grid);

    Vertex inverseRk3Integrate(Vertex newPosition, StaggeredVelocityGrid& grid);

protected:
    StaggeredVelocityGrid m_advectedVelocity;
    SdfGrid m_advectedSdf;
    Grid2d<float> m_uWeights;
    Grid2d<float> m_vWeights;
};

#endif // NBFLIPSOLVER_H
