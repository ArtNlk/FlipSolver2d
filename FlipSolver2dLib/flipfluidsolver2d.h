#ifndef FLIPSOLVER_H
#define FLIPSOLVER_H

#include "flipsolverbase.h"

class FlipFluidSolver : public FlipSolverBase
{
public:
    FlipFluidSolver(int extrapRadius = 1, bool vonNeumannNeighbors = false);

    ~FlipFluidSolver() override = default;

    void step() override;
};

#endif // FLIPSOLVER_H
