#ifndef MULTIFLIPSOLVER_H
#define MULTIFLIPSOLVER_H

#include "flipsolver2d.h"

class MultiflipSolver : public FlipSolver
{
public:
    MultiflipSolver(int extrapRadius = 1, bool vonNeumannNeighbors = false);

    void step() override;

    void updateSdf() override;

    void combineSdf();

    void extrapolateSdf(Grid2d<float>& sdfGrid);

    void countParticles() override;

    void reseedParticles() override;

    void seedInitialFluid() override;

    void updateMaterialsFromParticles() override;

    void bumpParticles();
};

#endif // MULTIFLIPSOLVER_H
