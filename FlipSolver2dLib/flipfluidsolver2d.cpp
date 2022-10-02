#include "flipfluidsolver2d.h"

FlipFluidSolver::FlipFluidSolver(int extrapRadius, bool vonNeumannNeighbors) :
    FlipSolver(extrapRadius, vonNeumannNeighbors)
{

}

void FlipFluidSolver::step()
{
    Grid2d<int> particleCounts(m_grid.sizeI(), m_grid.sizeJ());
    m_grid.updateLinearFluidViscosityMapping();
    countParticles(particleCounts);
    reseedParticles(particleCounts);
    updateMaterialsFromParticles(particleCounts);
    particleToGrid();
    extrapolateVelocityField(1);
    Grid2d<float> prevU = m_grid.velocityGridU();
    Grid2d<float> prevV = m_grid.velocityGridV();
    applyBodyForces();
    project();
    updateVelocityFromSolids();
    applyViscosity();
    project();
    extrapolateVelocityField(1);
    particleUpdate(prevU, prevV);
    advect();
}
