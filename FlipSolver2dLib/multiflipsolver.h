#ifndef MULTIFLIPSOLVER_H
#define MULTIFLIPSOLVER_H

#include "flipsolver2d.h"

class MultiflipSolver : public FlipSolver
{
public:
    MultiflipSolver(int extrapRadius = 1, bool vonNeumannNeighbors = false);

    void step() override;
};

#endif // MULTIFLIPSOLVER_H
