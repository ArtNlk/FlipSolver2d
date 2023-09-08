#ifndef FLIPFIRESOLVER_H
#define FLIPFIRESOLVER_H

#include "flipsmokesolver.h"

struct FireSolverParameters : SmokeSolverParameters
{
    float ignitionTemperature;
    float burnRate;
    float smokeProportion;
    float heatProportion;
    float divergenceProportion;
};

class FlipFireSolver : public FlipSmokeSolver
{
public:
    FlipFireSolver(const FireSolverParameters* p);

protected:
    void afterTransfer() override;
    void combustionUpdate();
    void combustionUpdateThread(Range range);
    void centeredParamsToGrid() override;
    void reseedParticles() override;
    void particleUpdate() override;

    void initAdditionalParameters() override;

    Grid2d<float> m_fuel;
    Grid2d<float> m_oxidizer;

    size_t m_fuelPropertyIndex;

    float m_ignitionTemperature;
    float m_burnRate;
    float m_smokeProportion;
    float m_heatProportion;
    float m_divergenceProportion;
};

#endif // FLIPFIRESOLVER_H
