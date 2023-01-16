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
    updateMaterials();
    particleToGrid();
    extrapolateVelocityField(m_grid.fluidVelocityGridU(),m_grid.knownFluidFlagsGridU(),10);
    extrapolateVelocityField(m_grid.fluidVelocityGridV(),m_grid.knownFluidFlagsGridV(),10);

    m_grid.savedFluidVelocityGrid().velocityGridU() = m_grid.fluidVelocityGridU();
    m_grid.savedFluidVelocityGrid().velocityGridV() = m_grid.fluidVelocityGridV();
    applyBodyForces();
    project();
    updateVelocityFromSolids();
    //applyViscosity();
    //project();
    extrapolateVelocityField(m_grid.fluidVelocityGridU(),m_grid.knownFluidFlagsGridU(),10);
    extrapolateVelocityField(m_grid.fluidVelocityGridV(),m_grid.knownFluidFlagsGridV(),10);

    particleUpdate();
    advect();
}
