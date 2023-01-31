#ifndef FLIPSMOKESOLVER_H
#define FLIPSMOKESOLVER_H

#include "dynamicuppertriangularsparsematrix.h"
#include "flipsolver2d.h"

class FlipSmokeSolver : public FlipSolver
{
public:
    FlipSmokeSolver(int sizeI, int sizeJ, int extrapRadius = 1, bool vonNeumannNeighbors = false);

    void applyBodyForces() override;

    void centeredParamsToGrid() override;

    void calcPressureRhs(std::vector<double> &rhs) override;

    void applyPressuresToVelocityField(std::vector<double> &pressures) override;

    DynamicUpperTriangularSparseMatrix getPressureProjectionMatrix() override;

    const Grid2d<float> smokeConcentration() const;

    const Grid2d<float> temperature() const;

protected:
    Grid2d<float> m_temperature;
    Grid2d<float> m_smokeConcentration;
};

#endif // FLIPSMOKESOLVER_H
