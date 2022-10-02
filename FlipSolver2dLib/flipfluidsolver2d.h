#ifndef FLIPFLUIDSOLVER_H
#define FLIPFLUIDSOLVER_H

#include "flipsolver2d.h"

class FlipFluidSolver : public FlipSolver
{
public:
    FlipFluidSolver(int extrapRadius = 1, bool vonNeumannNeighbors = false);

    ~FlipFluidSolver() override = default;

    void step() override;
};

#endif // FLIPFLUIDSOLVER_H
