#ifndef FLIPSMOKESOLVER_H
#define FLIPSMOKESOLVER_H

#include "dynamicuppertriangularsparsematrix.h"
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

    void calcPressureRhs(std::vector<double> &rhs) override;

    void particleUpdate() override;

    void afterTransfer() override;

    void reseedParticles() override;

    void applyPressuresToVelocityField(std::vector<double> &pressures) override;

    void step() override;

    void advect() override;

    LinearSolver::MatElementProvider getPressureMatrixElementProvider() override;

    LinearSolver::SparseMatRowElements getMatFreeElementForLinIdx(unsigned int i);

    DynamicUpperTriangularSparseMatrix getPressureProjectionMatrix() override;

protected:
    void eulerAdvectionThread(Range range, Vertex offset, const Grid2d<float>& inputGrid, Grid2d<float>& outputGrid);

    Vertex inverseRk4Integrate(Vertex newPosition, StaggeredVelocityGrid& grid);

    void combineAdvectedGrids();

    Grid2d<float> m_temperature;
    Grid2d<float> m_smokeConcentration;

    StaggeredVelocityGrid m_advectedVelocity;

    size_t m_temperatureIndex;
    size_t m_concentrationIndex;

    float m_ambientTemperature;
    float m_temperatureDecayRate;
    float m_concentrationDecayRate;
    float m_buoyancyFactor;
    float m_sootFactor;
};

#endif // FLIPSMOKESOLVER_H
