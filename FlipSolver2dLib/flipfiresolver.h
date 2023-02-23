#ifndef FLIPFIRESOLVER_H
#define FLIPFIRESOLVER_H

#include "flipsmokesolver.h"

class FlipFireSolver : public FlipSmokeSolver
{
public:
    FlipFireSolver(int sizeI, int sizeJ, int extrapRadius = 1, bool vonNeumannNeighbors = false);

protected:
    void afterTransfer() override;
    void combustionUpdate();
    void centeredParamsToGrid() override;
    void reseedParticles() override;
    void particleUpdate() override;

    Grid2d<float> m_fuel;
    Grid2d<float> m_oxidizer;
};

#endif // FLIPFIRESOLVER_H
