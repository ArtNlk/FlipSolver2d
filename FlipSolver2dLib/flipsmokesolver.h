#ifndef FLIPSMOKESOLVER_H
#define FLIPSMOKESOLVER_H

#include <Eigen/Dense>
#include "dynamicmatrix.h"
#include "flipsolver2d.h"

struct SmokeSolverParameters : FlipSolverParameters
{
    float ambientTemperature;
    float temperatureDecayRate;
    float concentrationDecayRate;
    float buoyancyFactor;
    float sootFactor;
};

class FlipSmokeSolver : public FlipSolver
{
public:
    FlipSmokeSolver(const SmokeSolverParameters* p);

    const Grid2d<float> smokeConcentration() const;

    const Grid2d<float> temperature() const;

    void initAdditionalParameters() override;

protected:

    void applyBodyForces() override;

    void centeredParamsToGrid() override;

    void calcPressureRhs(Eigen::VectorXd &rhs) override;

    void particleUpdate() override;

    void afterTransfer() override;

    void reseedParticles() override;

    void gridUpdate() override;

    void eulerAdvectParameters() override;

    void applyPressuresToVelocityField(Eigen::VectorXd &pressures) override;

    void seedInitialFluid() override;

    LinearSolver::MatElementProvider getPressureMatrixElementProvider() override;

    LinearSolver::SparseMatRowElements getMatFreeElementForLinIdx(unsigned int i);
    
    Eigen::SparseMatrix<double,Eigen::RowMajor> getPressureProjectionMatrix() override;

    void centeredParamsToGridThread(Range r, Grid2d<float>& cWeights);

    Grid2d<float> m_temperature;
    Grid2d<float> m_smokeConcentration;

    size_t m_temperatureIndex;
    size_t m_concentrationIndex;

    float m_ambientTemperature;
    float m_temperatureDecayRate;
    float m_concentrationDecayRate;
    float m_buoyancyFactor;
    float m_sootFactor;
};

#endif // FLIPSMOKESOLVER_H
