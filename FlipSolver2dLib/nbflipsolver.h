#ifndef NBFLIPSOLVER_H
#define NBFLIPSOLVER_H

#include "flipsolver2d.h"
#include "sdfgrid.h"
#include "staggeredvelocitygrid.h"

struct NBFlipParameters : FlipSolverParameters
{

};

class NBFlipSolver : public FlipSolver
{
public:
    NBFlipSolver(const NBFlipParameters* p);

protected:
    void step() override;

    void advect() override;

    void afterTransfer() override;

    void particleVelocityToGrid() override;

    void particleVelocityToGridThread(Range r, Grid2d<float>& uWeights, Grid2d<float>& vWeights);

    void reseedParticles() override;

    void firstFrameInit() override;

    void gridUpdate() override;

    void initialFluidSeed();

    void fluidSdfFromInitialFluid();

    void updateSdfFromSources();

    void combineAdvectedGrids();

    void combineLevelset();

    void combineVelocityGrid();
    
    Vertex inverseRk4Integrate(Vertex newPosition, StaggeredVelocityGrid& grid);

    const StaggeredVelocityGrid &advectedVelocityGrid() const;

protected:

    void centeredParamsToGrid() override;

    void centeredParamsToGridThread(Range r, Grid2d<float> &cWeights);

    StaggeredVelocityGrid m_advectedVelocity;
    SdfGrid m_advectedSdf;
    Grid2d<float> m_uWeights;
    Grid2d<float> m_vWeights;

    float m_narrowBand;
    float m_combinationBand;
    float m_resamplingBand;
};

#endif // NBFLIPSOLVER_H
