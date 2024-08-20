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

    void pruneNarrowBand();

    void initialFluidSeed();

    void fluidSdfFromInitialFluid();

    void updateGridFromSources();

    void combineAdvectedGrids();

    void combineLevelset();

    void combineVelocityGrid();

    void combineCenteredGrids();
    
    Vec3 inverseRk4Integrate(Vec3 newPosition, StaggeredVelocityGrid& grid);

    const StaggeredVelocityGrid &advectedVelocityGrid() const;

protected:

    inline float nbCombine(float grid, float particle, float sdf)
    {
        if(sdf > m_combinationBand)
        {
            return particle;
        }
        else
        {
            return grid;
        }
    };

    void centeredParamsToGrid() override;

    void centeredParamsToGridThread(Range r, Grid2d<float> &cWeights);

    Grid2d<float> m_advectedViscosity;
    StaggeredVelocityGrid m_advectedVelocity;
    SdfGrid m_advectedSdf;
    Grid2d<float> m_uWeights;
    Grid2d<float> m_vWeights;

    float m_narrowBand;
    float m_combinationBand;
    float m_resamplingBand;
};

#endif // NBFLIPSOLVER_H
