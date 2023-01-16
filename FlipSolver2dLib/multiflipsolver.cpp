#include "multiflipsolver.h"
#include "flipsolver2d.h"
#include "fluidcell.h"
#include "geometry2d.h"
#include "grid2d.h"
#include "logger.h"
#include "mathfuncs.h"
#include "simsettings.h"
#include "staggeredvelocitygrid.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

MultiflipSolver::MultiflipSolver(int extrapRadius, bool vonNeumannNeighbors) :
    FlipSolver(extrapRadius, vonNeumannNeighbors)
{

}

void MultiflipSolver::calcPressureRhs(std::vector<double> &rhs)
{
    double scale = 1.f/SimSettings::dx();
    Grid2d<float> curvatureGrid = simmath::calculateCenteredGridCurvature(m_grid.fluidSdfGrid());
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if (!m_grid.isSolid(i,j))
            {
                int linearIdx = m_grid.linearIndex(i,j);
                float weightedUCurrent = getWeightedVelocityUSample(i,j);
                float weightedVCurrent = getWeightedVelocityVSample(i,j);
                float weightedUIp1 = getWeightedVelocityUSample(i+1,j);
                float weightedVJp1 = getWeightedVelocityVSample(i,j+1);

                rhs[linearIdx] = -scale * (weightedUIp1-weightedUCurrent
                                        +weightedVJp1-weightedVCurrent);
                double surfaceTensionScale = SimSettings::stepDt() / (SimSettings::dx() * SimSettings::dx());
                double surfaceTensionFactor = SimSettings::surfaceTensionFactor() * curvatureGrid.at(i,j);

                if(m_grid.fluidSdfGrid().at(i,j) * m_grid.fluidSdfGrid().at(i-1,j) < 0.f)
                {
                    rhs[linearIdx] += surfaceTensionFactor * surfaceTensionScale * (1.f/getWeightedDensityForUSample(i,j));
                }
                if(m_grid.fluidSdfGrid().at(i,j) * m_grid.fluidSdfGrid().at(i,j-1) < 0.f)
                {
                    rhs[linearIdx] += surfaceTensionFactor * surfaceTensionScale * (1.f/getWeightedDensityForVSample(i,j));
                }

                if(m_grid.isSolid(i-1,j))
                {
                    rhs[linearIdx] -= scale * static_cast<double>(getWeightedVelocityUSample(i,j) - 0);
                }
                if(m_grid.isSolid(i+1,j))
                {
                    rhs[linearIdx] += scale * static_cast<double>(getWeightedVelocityUSample(i+1,j) - 0);
                }

                if(m_grid.isSolid(i,j-1))
                {
                    rhs[linearIdx] -= scale * static_cast<double>(getWeightedVelocityVSample(i,j) - 0);
                }
                if(m_grid.isSolid(i,j+1))
                {
                    rhs[linearIdx] += scale * static_cast<double>(getWeightedVelocityVSample(i,j+1) - 0);
                }

                m_grid.testGrid().at(i,j) = rhs[linearIdx];
                if(std::isnan(rhs[linearIdx]) || std::isinf(rhs[linearIdx]))
                {
                    std::cout << "Nan or inf in multiflip pressure RHS!\n";
                }
            }
        }
    }
}

DynamicUpperTriangularSparseMatrix MultiflipSolver::getPressureProjectionMatrix()
{
    DynamicUpperTriangularSparseMatrix output = DynamicUpperTriangularSparseMatrix(m_grid.cellCount(),7);

    //To be multiplied by 1/weighted density
    double preScale = SimSettings::stepDt() / (SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(!m_grid.isSolid(i,j))
            {
                double scaleU = preScale * (1.0 / getWeightedDensityForUSample(i,j));
                double scaleV = preScale * (1.0 / getWeightedDensityForVSample(i,j));
                double scaleUip1 = preScale * (1.0 / getWeightedDensityForUSample(i+1,j));
                double scaleVjp1 = preScale * (1.0 / getWeightedDensityForVSample(i,j+1));
                int rowIndex = m_grid.linearIndex(i,j);
                int im1Idx = m_grid.linearIndex(i-1,j);
                int jm1Idx = m_grid.linearIndex(i,j-1);

                if(std::isnan(scaleU) || std::isinf(scaleU)
                        ||std::isnan(scaleV) || std::isinf(scaleV)
                        ||std::isnan(scaleUip1) || std::isinf(scaleUip1)
                        ||std::isnan(scaleVjp1) || std::isinf(scaleVjp1))
                {
                    std::cout << "Nan in multiflip scales!" << '\n';
                }

                //X Neighbors
                if(!m_grid.isSolid(i-1,j) && m_grid.inBounds(i-1,j))
                {
                    output.addToAdiag(i,j,scaleU,  m_grid);
                }

                if(!m_grid.isSolid(i+1,j) && m_grid.inBounds(i+1,j))
                {
                    output.addToAdiag(i,j,scaleU,  m_grid);
                    output.setAx(i,j,-scaleUip1,  m_grid);
                }
//                if( m_grid.isFluid(i+1,j))
//                {
//                    output.addToAdiag(i,j,scaleU,  m_grid);
//                    output.setAx(i,j,-scaleUip1,  m_grid);
//                } else if( m_grid.isEmpty(i+1,j))
//                {
//                    output.addToAdiag(i,j,scaleU,  m_grid);
//                }

                //Y Neighbors
                if(!m_grid.isSolid(i,j-1) && m_grid.inBounds(i,j-1))
                {
                    output.addToAdiag(i,j,scaleV,  m_grid);
                }

                if(!m_grid.isSolid(i,j+1) && m_grid.inBounds(i,j+1))
                {
                    output.addToAdiag(i,j,scaleV,  m_grid);
                    output.setAy(i,j,-scaleVjp1,  m_grid);
                }
//                if( m_grid.isFluid(i,j+1))
//                {
//                    output.addToAdiag(i,j,scaleV,  m_grid);
//                    output.setAy(i,j,-scaleVjp1,  m_grid);
//                } else if( m_grid.isEmpty(i,j+1))
//                {
//                    output.addToAdiag(i,j,scaleV,  m_grid);
//                }
            }
        }
    }

    return output;
}

void MultiflipSolver::step()
{
    static int stepCount = 0;
    stepCount++;
    m_grid.updateLinearFluidViscosityMapping();
    countParticles();
    reseedParticles();
    updateSdf();
    combineSdf();
//    float maxSdf = *std::max_element(m_grid.fluidSdfGrid().data().begin(),m_grid.fluidSdfGrid().data().end());
//    if(maxSdf > 10000.f)
//    {
//        std::cout << "Bad sdf!";
//    }
    bumpParticles();
    updateMaterials();
    particleToGrid();
    extrapolateVelocityField(m_grid.fluidVelocityGridU(),m_grid.knownFluidFlagsGridU(),1e6);
    extrapolateVelocityField(m_grid.fluidVelocityGridV(),m_grid.knownFluidFlagsGridV(),1e6);
    extrapolateVelocityField(m_grid.airVelocityGridU(),m_grid.knownAirFlagsGridU(),1e6);
    extrapolateVelocityField(m_grid.airVelocityGridV(),m_grid.knownAirFlagsGridV(),1e6);

    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.isFluid(i,j))
            {
                m_grid.airVelocityGridU().at(i,j) = m_grid.fluidVelocityGridU().at(i,j);
                m_grid.airVelocityGridV().at(i,j) = m_grid.fluidVelocityGridV().at(i,j);
            }
        }
    }

    m_grid.savedFluidVelocityGrid().velocityGridU() = m_grid.fluidVelocityGridU();
    m_grid.savedFluidVelocityGrid().velocityGridV() = m_grid.fluidVelocityGridV();
    m_grid.savedAirVelocityGrid().velocityGridU() = m_grid.airVelocityGridU();
    m_grid.savedAirVelocityGrid().velocityGridV() = m_grid.airVelocityGridV();
    applyBodyForces();
    //return;
    project();
//    if(stepCount > 2)
//    {
//        return;
//    }
    //updateVelocityFromSolids();
    //applyViscosity();
    //project();
    //return;
    extrapolateVelocityField(m_grid.fluidVelocityGridU(),m_grid.knownFluidFlagsGridU(),1e6);
    extrapolateVelocityField(m_grid.fluidVelocityGridV(),m_grid.knownFluidFlagsGridV(),1e6);
    extrapolateVelocityField(m_grid.airVelocityGridU(),m_grid.knownAirFlagsGridU(),1e6);
    extrapolateVelocityField(m_grid.airVelocityGridV(),m_grid.knownAirFlagsGridV(),1e6);
    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.isFluid(i,j))
            {
                m_grid.airVelocityGridU().at(i,j) = m_grid.fluidVelocityGridU().at(i,j);
                m_grid.airVelocityGridV().at(i,j) = m_grid.fluidVelocityGridV().at(i,j);
            }
        }
    }
    particleUpdate();
    advect();
}

void MultiflipSolver::advect()
{
    //int maxSubsteps = -1;
    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
    {
        MarkerParticle &p = m_markerParticles[i];
        //int substepCount = 0;
        if(p.material == FluidMaterial::FLUID)
        {
            p.position = rk3Integrate(p.position,SimSettings::stepDt(), m_grid.fluidVelocityGrid());
        }
        else
        {
            p.position = rk3Integrate(p.position,SimSettings::stepDt(), m_grid.airVelocityGrid());
        }
        if(m_grid.solidSdfAt(p.position.x(),p.position.y()) < 0.f)
        {
            p.position = m_grid.closestSolidSurfacePoint(p.position);
        }
        //maxSubsteps = std::max(substepCount,maxSubsteps);
        int pI = simmath::integr(p.position.x());
        int pJ = simmath::integr(p.position.y());
        if(!m_grid.inBounds(pI,pJ) || m_grid.isSink(pI,pJ))
        {
            m_markerParticles.erase(markerParticles().begin() + i);
        }

    }
    //std::cout << "Multiflip advection done in max " << maxSubsteps << " substeps" << std::endl;
}

void MultiflipSolver::updateSdf()
{
    float particleRadius = SimSettings::particleScale()*SimSettings::dx();
    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            float fluidDistSqrd = std::numeric_limits<float>::max();
            float airDistSqrd = std::numeric_limits<float>::max();
            Vertex centerPoint = Vertex(static_cast<float>(i) + 0.5f,
                                        static_cast<float>(j) + 0.5f);
            for(MarkerParticle& p : m_markerParticles)
            {
                float diffX = p.position.x() - centerPoint.x();
                float diffY = p.position.y() - centerPoint.y();
                float newDistSqrd = diffX*diffX + diffY * diffY;
                if(p.material == FluidMaterial::FLUID)
                {
                    if(newDistSqrd < fluidDistSqrd)
                    {
                        fluidDistSqrd = newDistSqrd;
                    }
                }
                else
                {
                    if(newDistSqrd < airDistSqrd)
                    {
                        airDistSqrd = newDistSqrd;
                    }
                }
            }

            m_grid.fluidSdfGrid().at(i,j) = std::sqrt(fluidDistSqrd) * SimSettings::dx() - particleRadius;
            m_grid.airSdfGrid().at(i,j) = std::sqrt(airDistSqrd) * SimSettings::dx() - particleRadius;

        }
    }
}

void MultiflipSolver::combineSdf()
{
    extrapolateSdf(m_grid.airSdfGrid());
    extrapolateSdf(m_grid.fluidSdfGrid());
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            m_grid.fluidSdfGrid().at(i,j) = (m_grid.fluidSdfGrid().at(i,j)-m_grid.airSdfGrid().at(i,j))/2.f;
        }
    }
}

void MultiflipSolver::extrapolateSdf(Grid2d<float> &sdfGrid)
{
    Grid2d<float> normalDeriv(sdfGrid.sizeI(), sdfGrid.sizeJ(), 1e5f, OOBStrategy::OOB_EXTEND);
    Grid2d<bool> extrapFlags(sdfGrid.sizeI(), sdfGrid.sizeJ(), true, OOBStrategy::OOB_EXTEND);

    std::function<float (Grid2d<float> &, Vertex &, void*)> updateFunc = simmath::normalDerivLinearExapolationUpdate;

    for (int i = 0; i < sdfGrid.sizeI(); i++)
    {
        for (int j = 0; j < sdfGrid.sizeJ(); j++)
        {
            if(i == 0
            || i == (sdfGrid.sizeI() - 1)
            || j == 0
            || j == (sdfGrid.sizeJ() - 1))
            {
                extrapFlags.at(i,j) = false;
            }
            if(sdfGrid.at(i,j) >= 0.f)
            {
                Vertex grad = simmath::gradCenteredGrid(i,j,sdfGrid);
                Vertex normal = grad.normalized();
                normalDeriv.at(i,j) = normal.dot(grad);
                extrapFlags.at(i,j) = false;
            }
        }
    }
    simmath::fastSweep(normalDeriv,extrapFlags,updateFunc, nullptr);

    if(anyNanInf(normalDeriv.data()))
    {
        std::cout << "Nan/Inf in fast sweep normal derivative\n";
    }

    updateFunc = simmath::sdfLinearExapolationUpdate;

    simmath::fastSweep(sdfGrid,extrapFlags,updateFunc, &normalDeriv);

    if(anyNanInf(sdfGrid.data()))
    {
        std::cout << "Nan/Inf in fast sweep after extrapolation\n";
    }

    for (int i = 0; i < sdfGrid.sizeI(); i++)
    {
        for (int j = 0; j < sdfGrid.sizeJ(); j++)
        {
            if(extrapFlags.at(i,j))
            {
                sdfGrid.at(i,j) *= -1.f;
            }
        }
    }
}

void MultiflipSolver::particleToGrid()
{
    ASSERT(m_grid.fluidVelocityGridU().sizeI() == m_grid.fluidVelocityGridU().sizeI() && m_grid.fluidVelocityGridU().sizeJ() == m_grid.fluidVelocityGridU().sizeJ());
    ASSERT(m_grid.fluidVelocityGridV().sizeI() == m_grid.fluidVelocityGridV().sizeI() && m_grid.fluidVelocityGridV().sizeJ() == m_grid.fluidVelocityGridV().sizeJ());

    Grid2d<float> uFluidWeights(m_grid.fluidVelocityGridU().sizeI(),m_grid.fluidVelocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vFluidWeights(m_grid.fluidVelocityGridV().sizeI(),m_grid.fluidVelocityGridV().sizeJ(),1e-10f);
    Grid2d<float> uAirWeights(m_grid.airVelocityGridU().sizeI(),m_grid.airVelocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vAirWeights(m_grid.airVelocityGridV().sizeI(),m_grid.airVelocityGridV().sizeJ(),1e-10f);
    Grid2d<float> centeredWeights(m_grid.sizeI(),m_grid.sizeJ(),1e-10f);

    m_grid.fluidVelocityGridU().fill(0.f);
    m_grid.fluidVelocityGridV().fill(0.f);
    m_grid.airVelocityGridU().fill(0.f);
    m_grid.airVelocityGridV().fill(0.f);
    m_grid.viscosityGrid().fill(0.f);
    m_grid.knownFluidFlagsGridU().fill(false);
    m_grid.knownFluidFlagsGridV().fill(false);
    m_grid.knownAirFlagsGridU().fill(false);
    m_grid.knownAirFlagsGridV().fill(false);
    m_grid.knownFlagsCenteredParams().fill(false);
    m_grid.smokeConcentrationGrid().fill(0.f);
    m_grid.temperatureGrid().fill(0);

    for(MarkerParticle &p : m_markerParticles)
    {
        int i = simmath::integr(p.position.x());
        int j = simmath::integr(p.position.y());
        //Run over all cells that this particle might affect
        for (int iOffset = -3; iOffset < 3; iOffset++)
        {
            for (int jOffset = -3; jOffset < 3; jOffset++)
            {
                int iIdx = i+iOffset;
                int jIdx = j+jOffset;
                float weightU = simmath::quadraticBSpline(p.position.x() - static_cast<float>(iIdx),
                                                     p.position.y() - (static_cast<float>(jIdx) + 0.5f));

                float weightV = simmath::quadraticBSpline(p.position.x() - (static_cast<float>(iIdx) + 0.5f),
                                                     p.position.y() - static_cast<float>(jIdx));
                float weightCentered = simmath::quadraticBSpline(p.position.x() - static_cast<float>(iIdx),
                                                     p.position.y() - static_cast<float>(jIdx));
                if(uFluidWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    //if(std::abs(weightU) > 1e-9f)
                    {
                        if(p.material == FluidMaterial::FLUID)
                        {
                            uFluidWeights.at(iIdx,jIdx) += weightU;
                            m_grid.fluidVelocityGridU().at(iIdx,jIdx) +=
                                    weightU * (p.velocity.x() * SimSettings::dx());
                            m_grid.knownFluidFlagsGridU().at(iIdx,jIdx) = true;
                        }
                        else
                        {
                            uAirWeights.at(iIdx,jIdx) += weightU;
                            m_grid.airVelocityGridU().at(iIdx,jIdx) +=
                                    weightU * (p.velocity.x() * SimSettings::dx());
                            m_grid.knownAirFlagsGridU().at(iIdx,jIdx) = true;
                        }
                    }
                }

                if(vFluidWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    //if(std::abs(weightV) > 1e-9f)
                    {
                        if(p.material == FluidMaterial::FLUID)
                        {
                            vFluidWeights.at(iIdx,jIdx) += weightV;
                            m_grid.fluidVelocityGridV().at(iIdx,jIdx) +=
                                    weightV * (p.velocity.y() * SimSettings::dx());
                            m_grid.knownFluidFlagsGridV().at(iIdx,jIdx) = true;
                        }
                        else
                        {
                            vAirWeights.at(iIdx,jIdx) += weightV;
                            m_grid.airVelocityGridV().at(iIdx,jIdx) +=
                                    weightV * (p.velocity.y() * SimSettings::dx());
                            m_grid.knownAirFlagsGridV().at(iIdx,jIdx) = true;
                        }
                    }
                }

//                if(centeredWeights.inBounds(iIdx,jIdx))
//                {
//                    if(std::abs(weightCentered) > 1e-9f)
//                    {
//                        centeredWeights.at(iIdx,jIdx) += weightCentered;
//                        m_grid.viscosityGrid().at(iIdx,jIdx) += weightCentered * p.viscosity;
//                        m_grid.temperatureGrid().at(iIdx,jIdx) += weightCentered * p.temperature;
//                        m_grid.smokeConcentrationGrid().at(iIdx,jIdx) += weightCentered * p.smokeConcentrartion;
//                        m_grid.knownFlagsCenteredParams().at(iIdx,jIdx) = true;
//                    }
//                }
            }
        }
    }

    for (int i = 0; i < m_grid.sizeI() + 1; i++)
    {
        for (int j = 0; j < m_grid.sizeJ() + 1; j++)
        {
            if(m_grid.fluidVelocityGridU().inBounds(i,j))
            {
                m_grid.fluidVelocityGridU().at(i,j) /= uFluidWeights.at(i,j);
                m_grid.airVelocityGridU().at(i,j) /= uAirWeights.at(i,j);
            }
            if(m_grid.fluidVelocityGridV().inBounds(i,j))
            {
                m_grid.fluidVelocityGridV().at(i,j) /= vFluidWeights.at(i,j);
                m_grid.airVelocityGridV().at(i,j) /= vAirWeights.at(i,j);
            }
//            if(m_grid.knownFlagsCenteredParams().inBounds(i,j))
//            {
//                if(m_grid.knownFlagsCenteredParams().at(i,j))
//                {
//                    m_grid.viscosityGrid().at(i,j) /= centeredWeights.at(i,j);
//                    m_grid.smokeConcentrationGrid().at(i,j) /= centeredWeights.at(i,j);
//                    m_grid.temperatureGrid().at(i,j) /= centeredWeights.at(i,j);
//                }
//                else
//                {
//                    m_grid.temperatureGrid().at(i,j) = SimSettings::ambientTemp();
//                }
//            }
        }
    }
}

void MultiflipSolver::applyBodyForces()
{
    for (int i = 0; i < m_grid.sizeI() + 1; i++)
    {
        for (int j = 0; j < m_grid.sizeJ() + 1; j++)
        {
            if(m_grid.fluidVelocityGridU().inBounds(i,j))
            {
                m_grid.fluidVelocityGridU().at(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().x();
                m_grid.airVelocityGridU().at(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().x();
            }
            if(m_grid.fluidVelocityGridV().inBounds(i,j))
            {
                m_grid.fluidVelocityGridV().at(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().y();
                m_grid.airVelocityGridV().at(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().y();
            }
        }
    }
}

void MultiflipSolver::countParticles()
{
    m_grid.fluidParticleCountGrid().fill(0);
    m_grid.airParticleCountGrid().fill(0);
    for(MarkerParticle& p : m_markerParticles)
    {
        p.testValue = 0.f;
        int i = std::floor(p.position.x());
        int j = std::floor(p.position.y());
        if(p.material == FluidMaterial::FLUID)
        {
            m_grid.fluidParticleCountGrid().at(i,j) += 1;
        }
        else
        {
            m_grid.airParticleCountGrid().at(i,j) += 1;
        }
//        if(m_grid.isSolid(i,j))
//        {
//            std::cout << "Particle in solid at " << i << "," << j << '\n';
//            debug() << "Particle in solid at " << i << "," << j;
//        }
    }
}

void MultiflipSolver::reseedParticles()
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            int fluidParticleCount = m_grid.fluidParticleCountGrid().at(i,j);
            int airParticleCount = m_grid.airParticleCountGrid().at(i,j);
            if(fluidParticleCount > 20)
            {
                std::cout << "too many particles " << fluidParticleCount << " at " << i << ' ' << j;
            }
            //std::cout << particleCount << " at " << i << " , " << j << std::endl;
            int additionalFluidParticles = SimSettings::particlesPerCell() - fluidParticleCount;
            int additionalAirParticles = SimSettings::particlesPerCell() - airParticleCount;
            if(additionalFluidParticles <= 0)
            {
                continue;
            }
            if(m_grid.isSource(i,j))
            {
                for(int p = 0; p < additionalFluidParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_grid.fluidVelocityAt(pos);
                    int emitterId = m_grid.emitterId(i,j);
                    float viscosity = m_sources[emitterId].viscosity();
                    float conc = m_sources[emitterId].concentrartion();
                    float temp = m_sources[emitterId].temperature();
                    FluidMaterial mat = FluidMaterial::FLUID;
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc,mat});
                }
            }
            else if (m_grid.isEmpty(i,j))
            {
                for(int p = 0; p < additionalAirParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_grid.fluidVelocityAt(pos);
                    float viscosity = m_grid.viscosityAt(pos);
                    float conc = m_grid.smokeConcentrationAt(pos);
                    float temp = m_grid.temperatureAt(pos);
                    FluidMaterial mat = FluidMaterial::EMPTY;
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc,mat});
                }
            }
//            else if(m_grid.fluidSdfGrid().at(i,j) < SimSettings::particleScale() * SimSettings::dx())
//            {
//                for(int p = 0; p < additionalParticles; p++)
//                {
//                    Vertex pos = jitteredPosInCell(i,j);
//                    Vertex velocity = m_grid.fluidVelocityAt(pos);
//                    float viscosity = m_grid.viscosityAt(pos);
//                    float conc = m_grid.smokeConcentrationAt(pos);
//                    float temp = m_grid.temperatureAt(pos);
//                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc});
//                }
//            }
        }
    }
}

void MultiflipSolver::seedInitialFluid()
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.isStrictFluid(i,j) || m_grid.isEmpty(i,j))
            {
                FluidMaterial particleMaterial = m_grid.getMaterial(i,j);
                for(int p = 0; p < SimSettings::particlesPerCell(); p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_grid.fluidVelocityAt(pos);
                    float viscosity = m_grid.viscosityAt(pos);
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,0,0,particleMaterial});
                }
            }
        }
    }
}

void MultiflipSolver::updateMaterials()
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.fluidParticleCountGrid().at(i,j) != 0)
            {
                if(m_grid.isEmpty(i,j))
                {
                    m_grid.setMaterial(i,j,FluidMaterial::FLUID);
                }
            }
            else
            {
                FluidMaterial m = m_grid.getMaterial(i,j);
                if(m == FluidMaterial::FLUID)
                {
                    m_grid.setMaterial(i,j,FluidMaterial::EMPTY);
                }
            }
//            if(!m_grid.isSolid(i,j))
//            {
//                if(m_grid.fluidSdfGrid().at(i,j) < 0.f)
//                {
//                    m_grid.setMaterial(i,j,FluidMaterial::FLUID);
//                }
//                else
//                {
//                    m_grid.setMaterial(i,j,FluidMaterial::EMPTY);
//                }
//            }
        }
    }
}

void MultiflipSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double preScale = SimSettings::stepDt() / (SimSettings::dx());

    Grid2d<float> curvatureGrid = simmath::calculateCenteredGridCurvature(m_grid.fluidSdfGrid());

    float minFluidFraction = 0.1;
    float maxFluidFraction = 0.9;

    for (int i = m_grid.sizeI() - 1; i >= 0; i--)
    {
        for (int j = m_grid.sizeJ() - 1; j >= 0; j--)
        {
            int currentIdx = m_grid.linearIndex(i,j);
            int im1Idx = m_grid.linearIndex(i-1,j);
            int jm1Idx = m_grid.linearIndex(i,j-1);
            double surfaceTensionFactor = SimSettings::surfaceTensionFactor() * curvatureGrid.at(i,j);
            //m_grid.testGrid().at(i,j) = pressures[currentIdx]/100.f;
            //U part
            //m_grid.testGrid().at(i,j) = getWeightedDensityForUSample(i,j);
            if(!m_grid.isSolid(i-1,j) || !m_grid.isSolid(i,j))
            {
                if(m_grid.isSolid(i-1,j) || m_grid.isSolid(i,j))
                {
                    m_grid.setFluidU(i,j,0);//Solids are stationary
                    m_grid.setAirU(i,j,0);
                    m_grid.knownFluidFlagsGridU().at(i,j) = true;
                    m_grid.knownAirFlagsGridU().at(i,j) = true;
                }
                else
                {
                    float fluidFaceFraction = m_grid.getFaceFractionUSample(i,j);
                    float invWeightedDensity = 1.f/getWeightedDensityForUSample(i,j);
                    float pressureGrad = pressures[currentIdx]
                                        - (im1Idx != -1? pressures[im1Idx] : 0.f);
                    //m_grid.testGrid().at(i,j) = pressures[currentIdx];
                    if(fluidFaceFraction > minFluidFraction)
                    {
                        if(fluidFaceFraction < maxFluidFraction)
                        {//Update both
                            m_grid.fluidVelocityGrid().u(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            m_grid.airVelocityGrid().u(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            if(m_grid.fluidSdfGrid().at(i,j) * m_grid.fluidSdfGrid().at(i-1,j) < 0.f)
                            {
                                m_grid.fluidVelocityGrid().u(i,j) += preScale * invWeightedDensity * surfaceTensionFactor;
                                m_grid.airVelocityGrid().u(i,j) += preScale * invWeightedDensity * surfaceTensionFactor;
                            }
                            m_grid.knownFluidFlagsGridU().at(i,j) = true;
                            m_grid.knownAirFlagsGridU().at(i,j) = true;
                        }
                        else
                        {//Only fluid
                            m_grid.fluidVelocityGrid().u(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            if(m_grid.fluidSdfGrid().at(i,j) * m_grid.fluidSdfGrid().at(i-1,j) < 0.f)
                            {
                                m_grid.fluidVelocityGrid().u(i,j) += preScale * invWeightedDensity * surfaceTensionFactor;
                            }
                            m_grid.knownFluidFlagsGridU().at(i,j) = true;
                            m_grid.knownAirFlagsGridU().at(i,j) = false;
                        }
                    }
                    else
                    {//Only air
                        m_grid.airVelocityGrid().u(i,j) -= preScale * invWeightedDensity * pressureGrad;
                        if(m_grid.fluidSdfGrid().at(i,j) * m_grid.fluidSdfGrid().at(i-1,j) < 0.f)
                        {
                            m_grid.airVelocityGrid().u(i,j) += preScale * invWeightedDensity * surfaceTensionFactor;
                        }
                        m_grid.knownFluidFlagsGridU().at(i,j) = false;
                        m_grid.knownAirFlagsGridU().at(i,j) = true;
                    }
                }
            }
            else
            {
                m_grid.knownFluidFlagsGridU().at(i,j) = false;
                m_grid.knownAirFlagsGridU().at(i,j) = false;
            }

            //V part
            if(!m_grid.isSolid(i,j-1) || !m_grid.isSolid(i,j))
            {
                if(m_grid.isSolid(i,j-1) || m_grid.isSolid(i,j))
                {
                    m_grid.setFluidV(i,j,0);//Solids are stationary
                    m_grid.setAirV(i,j,0);
                    m_grid.knownFluidFlagsGridV().at(i,j) = true;
                    m_grid.knownAirFlagsGridV().at(i,j) = true;
                }
                else
                {
                    float fluidFaceFraction = m_grid.getFaceFractionVSample(i,j);
                    float invWeightedDensity = 1.f/getWeightedDensityForVSample(i,j);
                    float pressureGrad = pressures[currentIdx]
                                        - (jm1Idx != -1? pressures[jm1Idx] : 0.f);
                    if(fluidFaceFraction > minFluidFraction)
                    {
                        if(fluidFaceFraction < maxFluidFraction)
                        {//Update both
                            m_grid.fluidVelocityGrid().v(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            m_grid.airVelocityGrid().v(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            if(m_grid.fluidSdfGrid().at(i,j) * m_grid.fluidSdfGrid().at(i,j-1) < 0.f)
                            {
                                m_grid.fluidVelocityGrid().v(i,j) += preScale * invWeightedDensity * surfaceTensionFactor;
                                m_grid.airVelocityGrid().v(i,j) += preScale * invWeightedDensity * surfaceTensionFactor;
                            }
                            m_grid.knownFluidFlagsGridV().at(i,j) = true;
                            m_grid.knownAirFlagsGridV().at(i,j) = true;
                        }
                        else
                        {//Only fluid
                            m_grid.fluidVelocityGrid().v(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            if(m_grid.fluidSdfGrid().at(i,j) * m_grid.fluidSdfGrid().at(i,j-1) < 0.f)
                            {
                                m_grid.fluidVelocityGrid().v(i,j) += preScale * invWeightedDensity * surfaceTensionFactor;
                            }
                            m_grid.knownFluidFlagsGridV().at(i,j) = true;
                            m_grid.knownAirFlagsGridV().at(i,j) = false;
                        }
                    }
                    else
                    {//Only air
                        m_grid.airVelocityGrid().v(i,j) -= preScale * invWeightedDensity * pressureGrad;
                        if(m_grid.fluidSdfGrid().at(i,j) * m_grid.fluidSdfGrid().at(i,j-1) < 0.f)
                        {
                            m_grid.airVelocityGrid().v(i,j) += preScale * invWeightedDensity * surfaceTensionFactor;
                        }
                        m_grid.knownFluidFlagsGridV().at(i,j) = false;
                        m_grid.knownAirFlagsGridV().at(i,j) = true;
                    }
                }
            }
            else
            {
                m_grid.knownFluidFlagsGridV().at(i,j) = false;
                m_grid.knownAirFlagsGridV().at(i,j) = false;
            }
        }
    }

    if(anyNanInf(m_grid.fluidVelocityGridU().data()))
    {
        std::cout << "NaN or inf in U fluid!\n" << std::flush;
    }

    if(anyNanInf(m_grid.fluidVelocityGridV().data()))
    {
        std::cout << "NaN or inf in V fluid!\n" << std::flush;
    }

    if(anyNanInf(m_grid.airVelocityGridU().data()))
    {
        std::cout << "NaN or inf in U air!\n" << std::flush;
    }

    if(anyNanInf(m_grid.airVelocityGridV().data()))
    {
        std::cout << "NaN or inf in V air!\n" << std::flush;
    }
}

void MultiflipSolver::particleUpdate()
{
    Grid2d<float>& prevFluidU = m_grid.savedFluidVelocityGrid().velocityGridU();
    Grid2d<float>& prevFluidV = m_grid.savedFluidVelocityGrid().velocityGridV();
    Grid2d<float>& prevAirU = m_grid.savedAirVelocityGrid().velocityGridU();
    Grid2d<float>& prevAirV = m_grid.savedAirVelocityGrid().velocityGridV();

    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
    {
        MarkerParticle &p = m_markerParticles[i];
        if(p.material == FluidMaterial::FLUID)
        {
            Vertex oldFluidVelocity(simmath::lerpUGrid(p.position.x(),p.position.y(),prevFluidU)
                                    / SimSettings::dx(),
                               simmath::lerpVGrid(p.position.x(),p.position.y(),prevFluidV)
                                    / SimSettings::dx());
            Vertex newFluidVelocity = m_grid.fluidVelocityAt(p.position) / SimSettings::dx();
            p.velocity = SimSettings::picRatio() * newFluidVelocity +
                    (1.f-SimSettings::picRatio()) * (p.velocity + newFluidVelocity - oldFluidVelocity);
        }
        else
        {
            Vertex oldAirVelocity(simmath::lerpUGrid(p.position.x(),p.position.y(),prevAirU)
                                  / SimSettings::dx(),
                               simmath::lerpVGrid(p.position.x(),p.position.y(),prevAirV)
                                  / SimSettings::dx());
            Vertex newAirVelocity = m_grid.airVelocityAt(p.position) / SimSettings::dx();
            p.velocity = SimSettings::picRatio() * newAirVelocity +
                    (1.f-SimSettings::picRatio()) * (p.velocity + newAirVelocity - oldAirVelocity);
        }
    }
}

void MultiflipSolver::bumpParticles()
{
    const float escapeRadius = 1.5f * SimSettings::dx();
    const float particleRadius = SimSettings::particleScale() * SimSettings::dx();
    Grid2d<float> curvatureGrid = simmath::calculateCenteredGridCurvature(m_grid.fluidSdfGrid());
    for(MarkerParticle& p : m_markerParticles)
    {
        float sign = p.material == FluidMaterial::FLUID ? -1.f : 1.f;
        float sdfAtParticle = simmath::lerpCenteredGrid(p.position.x(),p.position.y(),m_grid.fluidSdfGrid());
        Vertex sdfGradAtParticle = simmath::gradCenteredGrid(p.position.x(),p.position.y(),m_grid.fluidSdfGrid());
        float d = sdfAtParticle / sdfGradAtParticle.distFromZero();
        if(sign*sdfAtParticle > -particleRadius && sign*sdfAtParticle < escapeRadius)
        {
            float curvature = simmath::lerpCenteredGrid(p.position, curvatureGrid);
            const float curvatureMin = 1.f/(SimSettings::dx()*4);
            const float curvatureMax = 1.f/(SimSettings::dx()*2);
            float negativeSignCurvature = -sign*curvature;
            float targetDist = 0.f;


            if(negativeSignCurvature > curvatureMin &&
                    negativeSignCurvature < curvatureMax)
            {
                targetDist = 0.f;
            }
            else if(negativeSignCurvature >= curvatureMax)
            {
                continue;
            }
            else
            {
                targetDist = particleRadius;
            }

            Vertex vectorToClosestSurfacePoint = p.position - m_grid.closesFluidSurfacePoint(p.position);
            float distanceToSurfacePoint = vectorToClosestSurfacePoint.distFromZero();

            if(sign * d < 0.f)
            {
                p.testValue = 1.f;
                p.position = p.position
                        -((targetDist + distanceToSurfacePoint)
                        *(vectorToClosestSurfacePoint/distanceToSurfacePoint))/SimSettings::dx();
            }
            else if (sign * d >= 0.f && sign * d < targetDist)
            {
                p.position = p.position
                        +((targetDist - distanceToSurfacePoint)
                        *(vectorToClosestSurfacePoint/distanceToSurfacePoint))/SimSettings::dx();
            }
        }
    }
}

double MultiflipSolver::getWeightedDensityForUSample(int i, int j)
{
    float sdfCurrent = m_grid.fluidSdfGrid().at(i,j) / SimSettings::dx();
    float sdfAtIm1 = m_grid.fluidSdfGrid().at(i-1,j) / SimSettings::dx();
    Vertex sdfGradAtSample = simmath::gradCenteredGrid(static_cast<float>(i) + 0.5f,
                                                       static_cast<float>(j) + 0.5f,
                                                       m_grid.fluidSdfGrid());
    float normDistCurr = sdfCurrent / sdfGradAtSample.distFromZero();

    Vertex sdfGradAtIm1 = simmath::gradCenteredGrid(static_cast<float>(i) - 0.5f,
                                                       static_cast<float>(j) + 0.5f,
                                                       m_grid.fluidSdfGrid());
    float normDistIm1 = sdfCurrent / sdfGradAtIm1.distFromZero();
    if(std::signbit(sdfCurrent) == std::signbit(sdfAtIm1))
    {
        if(sdfCurrent >= 0.f)
        {
            return SimSettings::airDensity();
        }
        else
        {
            return SimSettings::fluidDensity();
        }
    }

    float fluidFrac = 0.f;
    if(sdfCurrent < 0.f)
    {
        fluidFrac = -sdfCurrent;
    }
    else
    {
        fluidFrac = -sdfAtIm1;
    }
    //fluidFrac = m_grid.getFaceFractionUSample(i,j);
    fluidFrac = std::clamp(fluidFrac,0.f,1.f);
    float result = SimSettings::fluidDensity() * fluidFrac + (1.f - fluidFrac) * SimSettings::airDensity();
    return result;
}

double MultiflipSolver::getWeightedDensityForVSample(int i, int j)
{
    float sdfCurrent = m_grid.fluidSdfGrid().at(i,j) / SimSettings::dx();
    float sdfAtJm1 = m_grid.fluidSdfGrid().at(i,j-1) / SimSettings::dx();
    Vertex sdfGradAtSample = simmath::gradCenteredGrid(static_cast<float>(i) + 0.5f,
                                                       static_cast<float>(j) + 0.5f,
                                                       m_grid.fluidSdfGrid());
    float normDistCurr = sdfCurrent / sdfGradAtSample.distFromZero();

    Vertex sdfGradAtJm1 = simmath::gradCenteredGrid(static_cast<float>(i) + 0.5f,
                                                       static_cast<float>(j) - 0.5f,
                                                       m_grid.fluidSdfGrid());
    float normDistJm1 = sdfCurrent / sdfGradAtJm1.distFromZero();
    if(std::signbit(sdfCurrent) == std::signbit(sdfAtJm1))
    {
        if(sdfCurrent >= 0.f)
        {
            return SimSettings::airDensity();
        }
        else
        {
            return SimSettings::fluidDensity();
        }
    }

    float fluidFrac = 0.f;
    if(sdfCurrent < 0.f)
    {
        fluidFrac = -sdfCurrent;
    }
    else
    {
        fluidFrac = -sdfAtJm1;
    }
    //fluidFrac = m_grid.getFaceFractionVSample(i,j);
    fluidFrac = std::clamp(fluidFrac,0.f,1.f);
    float result = SimSettings::fluidDensity() * fluidFrac + (1.f - fluidFrac) * SimSettings::airDensity();
    return result;
}

double MultiflipSolver::getWeightedVelocityUSample(int i, int j)
{
    double faceFraction = m_grid.getFaceFractionUSample(i,j);
    double fluidVelocity = m_grid.getFluidU(i,j);
    double airVelocity = m_grid.getAirU(i,j);

    return fluidVelocity * faceFraction + airVelocity * (1.f - faceFraction);
}

double MultiflipSolver::getWeightedVelocityVSample(int i, int j)
{
    double faceFraction = m_grid.getFaceFractionVSample(i,j);
    double fluidVelocity = m_grid.getFluidV(i,j);
    double airVelocity = m_grid.getAirV(i,j);

    return fluidVelocity * faceFraction + airVelocity * (1.f - faceFraction);
}
