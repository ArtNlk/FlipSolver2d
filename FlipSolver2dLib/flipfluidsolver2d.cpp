#include "flipfluidsolver2d.h"

FlipFluidSolver::FlipFluidSolver(int extrapRadius, bool vonNeumannNeighbors) :
    FlipSolver(extrapRadius, vonNeumannNeighbors)
{

}

void FlipFluidSolver::step()
{
    updateSdf();
    m_grid.updateLinearFluidViscosityMapping();
    countParticles();
    reseedParticles();
    updateMaterialsFromParticles();
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
