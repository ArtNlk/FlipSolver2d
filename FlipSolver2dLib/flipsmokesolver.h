#ifndef FLIPSMOKESOLVER_H
#define FLIPSMOKESOLVER_H

#include "dynamicuppertriangularsparsematrix.h"
#include "flipsolver2d.h"

struct SmokeSolverParameters : FlipSolverParameters
{
    float ambientTemperature;
    float temperatureDecayRate;
    float concentrationDecayRate;
};

class FlipSmokeSolver : public FlipSolver
{
public:
    FlipSmokeSolver(const SmokeSolverParameters* p);

    const Grid2d<float> smokeConcentration() const;

    const Grid2d<float> temperature() const;

protected:
    void applyBodyForces() override;

    void centeredParamsToGrid() override;

    void calcPressureRhs(std::vector<double> &rhs) override;

    void particleUpdate() override;

    void afterTransfer() override;

    void reseedParticles() override;

    void applyPressuresToVelocityField(std::vector<double> &pressures) override;

    DynamicUpperTriangularSparseMatrix getPressureProjectionMatrix() override;

protected:
    Grid2d<float> m_temperature;
    Grid2d<float> m_smokeConcentration;

    float m_ambientTemperature;
    float m_temperatureDecayRate;
    float m_concentrationDecayRate;
};

#endif // FLIPSMOKESOLVER_H
