#ifndef FLIPSMOKESOLVER_H
#define FLIPSMOKESOLVER_H

#include "dynamicuppertriangularsparsematrix.h"
#include "flipsolver2d.h"

class FlipSmokeSolver : public FlipSolver
{
public:
    FlipSmokeSolver(int extrapRadius = 1, bool vonNeumannNeighbors = false);

    void step() override;

    void applyBodyForces() override;

    void particleToGrid() override;

    void project() override;

    void calcPressureRhs(std::vector<double> &rhs) override;

    void applyPressuresToVelocityField(std::vector<double> &pressures) override;

    DynamicUpperTriangularSparseMatrix getPressureProjectionMatrix() override;
};

#endif // FLIPSMOKESOLVER_H
