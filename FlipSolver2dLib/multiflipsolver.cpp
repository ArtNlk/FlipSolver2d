#include "multiflipsolver.h"
#include "flipsolver2d.h"
#include "materialgrid.h"
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

MultiflipSolver::MultiflipSolver(int sizeI, int sizeJ, int extrapRadius, bool vonNeumannNeighbors) :
    FlipSolver(sizeI, sizeJ, extrapRadius, vonNeumannNeighbors),
    m_airVelocityGrid(sizeI,sizeJ),
    m_savedAirVelocityGrid(sizeI, sizeJ),
    m_airSdf(sizeI, sizeJ),
    m_airParticleCounts(sizeI, sizeJ)
{

}

void MultiflipSolver::calcPressureRhs(std::vector<double> &rhs)
{
    double scale = 1.f/SimSettings::dx();
    Grid2d<float> curvatureGrid = simmath::calculateCenteredGridCurvature(m_fluidSdf);
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if (!m_materialGrid.isSolid(i,j))
            {
                int linearIdx = linearIndex(i,j);

                rhs[linearIdx] = -scale * divergenceAt(i,j);
                double surfaceTensionScale = SimSettings::stepDt() / (SimSettings::dx() * SimSettings::dx());
                Vertex uSamplePos(0.5f + i, j);
                Vertex vSamplePos(i, 0.5f + j);
                double surfaceTensionFactorU = SimSettings::surfaceTensionFactor()
                                                * curvatureGrid.interpolateAt(uSamplePos);
                double surfaceTensionFactorV = SimSettings::surfaceTensionFactor()
                                                * curvatureGrid.interpolateAt(vSamplePos);

                if(m_fluidSdf.getAt(i,j) * m_fluidSdf.getAt(i-1,j) < 0.f ||
                    m_fluidSdf.getAt(i,j) * m_fluidSdf.getAt(i,j-1) < 0.f)
                {
                    rhs[linearIdx] += surfaceTensionFactorU * surfaceTensionScale * (1.f/getWeightedDensityForUSample(i,j));
                    //m_testGrid.at(i,j) = 1.f;
                    rhs[linearIdx] += surfaceTensionFactorV * surfaceTensionScale * (1.f/getWeightedDensityForVSample(i,j));
                    //m_testGrid.at(i,j) = 1.f;
                }

                if(m_materialGrid.isSolid(i-1,j))
                {
                    rhs[linearIdx] -= scale * static_cast<double>(getWeightedVelocityUSample(i,j) - 0.);
                }
                if(m_materialGrid.isSolid(i+1,j))
                {
                    rhs[linearIdx] += scale * static_cast<double>(getWeightedVelocityUSample(i+1,j) - 0.);
                }

                if(m_materialGrid.isSolid(i,j-1))
                {
                    rhs[linearIdx] -= scale * static_cast<double>(getWeightedVelocityVSample(i,j) - 0.);
                }
                if(m_materialGrid.isSolid(i,j+1))
                {
                    rhs[linearIdx] += scale * static_cast<double>(getWeightedVelocityVSample(i,j+1) - 0.);
                }

                //m_testGrid.at(i,j) = rhs[linearIdx];
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
    DynamicUpperTriangularSparseMatrix output = DynamicUpperTriangularSparseMatrix(cellCount(),7);

    //To be multiplied by 1/weighted density
    double preScale = SimSettings::stepDt() / (SimSettings::dx() * SimSettings::dx());

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(!m_materialGrid.isSolid(i,j))
            {
                double scaleU = preScale * (1.0 / getWeightedDensityForUSample(i,j));
                double scaleV = preScale * (1.0 / getWeightedDensityForVSample(i,j));
                double scaleUip1 = preScale * (1.0 / getWeightedDensityForUSample(i+1,j));
                double scaleVjp1 = preScale * (1.0 / getWeightedDensityForVSample(i,j+1));
                int rowIndex = linearIndex(i,j);
                int im1Idx = linearIndex(i-1,j);
                int jm1Idx = linearIndex(i,j-1);

                if(std::isnan(scaleU) || std::isinf(scaleU)
                        ||std::isnan(scaleV) || std::isinf(scaleV)
                        ||std::isnan(scaleUip1) || std::isinf(scaleUip1)
                        ||std::isnan(scaleVjp1) || std::isinf(scaleVjp1))
                {
                    std::cout << "Nan in multiflip scales!" << '\n';
                }

                //X Neighbors
                if(!m_materialGrid.isSolid(i-1,j) && inBounds(i-1,j))
                {
                    output.addToAdiag(i,j,scaleU,  indexer);
                }

                if(!m_materialGrid.isSolid(i+1,j) && inBounds(i+1,j))
                {
                    output.addToAdiag(i,j,scaleU,  indexer);
                    output.setAx(i,j,-scaleUip1,  indexer);
                }
//                if( m_materialGrid.isFluid(i+1,j))
//                {
//                    output.addToAdiag(i,j,scaleU,  indexer);
//                    output.setAx(i,j,-scaleUip1,  indexer);
//                } else if( m_materialGrid.isEmpty(i+1,j))
//                {
//                    output.addToAdiag(i,j,scaleU,  indexer);
//                }

                //Y Neighbors
                if(!m_materialGrid.isSolid(i,j-1) && inBounds(i,j-1))
                {
                    output.addToAdiag(i,j,scaleV,  indexer);
                }

                if(!m_materialGrid.isSolid(i,j+1) && inBounds(i,j+1))
                {
                    output.addToAdiag(i,j,scaleV,  indexer);
                    output.setAy(i,j,-scaleVjp1,  indexer);
                }
//                if( m_materialGrid.isFluid(i,j+1))
//                {
//                    output.addToAdiag(i,j,scaleV,  indexer);
//                    output.setAy(i,j,-scaleVjp1,  indexer);
//                } else if( m_materialGrid.isEmpty(i,j+1))
//                {
//                    output.addToAdiag(i,j,scaleV,  indexer);
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
    updateLinearFluidViscosityMapping();
    countParticles();
    reseedParticles();
    updateSdf();
    combineSdf();
//    float maxSdf = *std::max_element(m_fluidSdf.data().begin(),m_fluidSdf.data().end());
//    if(maxSdf > 10000.f)
//    {
//        std::cout << "Bad sdf!";
//    }
    bumpParticles();
    updateMaterials();
    particleToGrid();
    m_fluidVelocityGrid.extrapolate(100);
    m_airVelocityGrid.extrapolate(100);

//    for(int i = 0; i < m_sizeI; i++)
//    {
//        for(int j = 0; j < m_sizeJ; j++)
//        {
//            if(m_materialGrid.isFluid(i,j))
//            {
//                m_airVelocityGrid.u(i,j) = m_fluidVelocityGrid.u(i,j);
//                m_airVelocityGrid.v(i,j) = m_fluidVelocityGrid.v(i,j);
//            }
//        }
//    }

    m_savedFluidVelocityGrid = m_fluidVelocityGrid;
    m_savedAirVelocityGrid = m_airVelocityGrid;
    applyBodyForces();
    //return;
    std::vector<double> rhs(cellCount(),0.0);
    project();
    calcPressureRhs(rhs);
//    if(stepCount > 2)
//    {
//        return;
//    }
    //updateVelocityFromSolids();
    //applyViscosity();
    //project();
    //return;
    m_fluidVelocityGrid.extrapolate(10);
    m_airVelocityGrid.extrapolate(10);
//    for(int i = 0; i < m_sizeI; i++)
//    {
//        for(int j = 0; j < m_sizeJ; j++)
//        {
//            if(m_materialGrid.isFluid(i,j))
//            {
//                m_airVelocityGrid.u(i,j) = m_fluidVelocityGrid.u(i,j);
//                m_airVelocityGrid.v(i,j) = m_fluidVelocityGrid.v(i,j);
//            }
//        }
//    }
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
            p.position = rk3Integrate(p.position,SimSettings::stepDt(), m_fluidVelocityGrid);
        }
        else
        {
            p.position = rk3Integrate(p.position,SimSettings::stepDt(), m_airVelocityGrid);
        }
        if(m_solidSdf.interpolateAt(p.position.x(),p.position.y()) < 0.f)
        {
            p.position = m_solidSdf.closestSurfacePoint(p.position);
        }
        //maxSubsteps = std::max(substepCount,maxSubsteps);
        int pI = simmath::integr(p.position.x());
        int pJ = simmath::integr(p.position.y());
        if(!inBounds(pI,pJ) || m_materialGrid.isSink(pI,pJ))
        {
            m_markerParticles.erase(markerParticles().begin() + i);
        }

    }
    //std::cout << "Multiflip advection done in max " << maxSubsteps << " substeps" << std::endl;
}

void MultiflipSolver::updateSdf()
{
    float particleRadius = SimSettings::particleScale()*SimSettings::dx();
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
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

            m_fluidSdf.at(i,j) = std::sqrt(fluidDistSqrd) * SimSettings::dx() - particleRadius;
            m_airSdf.at(i,j) = std::sqrt(airDistSqrd) * SimSettings::dx() - particleRadius;
        }
    }
}

void MultiflipSolver::combineSdf()
{
    extrapolateSdf(m_airSdf);
    extrapolateSdf(m_fluidSdf);
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            m_fluidSdf.at(i,j) = (m_fluidSdf.at(i,j)-m_airSdf.at(i,j))/2.f;
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

void MultiflipSolver::particleVelocityToGrid()
{
    Grid2d<float> uFluidWeights(m_fluidVelocityGrid.velocityGridU().sizeI(),m_fluidVelocityGrid.velocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vFluidWeights(m_fluidVelocityGrid.velocityGridV().sizeI(),m_fluidVelocityGrid.velocityGridV().sizeJ(),1e-10f);
    Grid2d<float> uAirWeights(m_airVelocityGrid.velocityGridU().sizeI(),m_airVelocityGrid.velocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vAirWeights(m_airVelocityGrid.velocityGridV().sizeI(),m_airVelocityGrid.velocityGridV().sizeJ(),1e-10f);

    m_fluidVelocityGrid.velocityGridU().fill(0.f);
    m_fluidVelocityGrid.velocityGridV().fill(0.f);
    m_airVelocityGrid.velocityGridU().fill(0.f);
    m_airVelocityGrid.velocityGridV().fill(0.f);

    m_fluidVelocityGrid.uSampleValidityGrid().fill(false);
    m_fluidVelocityGrid.vSampleValidityGrid().fill(false);
    m_airVelocityGrid.uSampleValidityGrid().fill(false);
    m_airVelocityGrid.vSampleValidityGrid().fill(false);

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
                if(uFluidWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    //if(std::abs(weightU) > 1e-9f)
                    {
                        if(p.material == FluidMaterial::FLUID)
                        {
                            uFluidWeights.at(iIdx,jIdx) += weightU;
                            m_fluidVelocityGrid.u(iIdx,jIdx) +=
                                    weightU * (p.velocity.x() * SimSettings::dx());
                            m_fluidVelocityGrid.setUValidity(iIdx,jIdx,true);
                        }
                        else
                        {
                            uAirWeights.at(iIdx,jIdx) += weightU;
                            m_airVelocityGrid.u(iIdx,jIdx) +=
                                    weightU * (p.velocity.x() * SimSettings::dx());
                            m_airVelocityGrid.setUValidity(iIdx,jIdx,true);
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
                            m_fluidVelocityGrid.v(iIdx,jIdx) +=
                                    weightV * (p.velocity.y() * SimSettings::dx());
                            m_fluidVelocityGrid.setVValidity(iIdx,jIdx,true);
                        }
                        else
                        {
                            vAirWeights.at(iIdx,jIdx) += weightV;
                            m_airVelocityGrid.v(iIdx,jIdx) +=
                                    weightV * (p.velocity.y() * SimSettings::dx());
                            m_airVelocityGrid.setVValidity(iIdx,jIdx,true);
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(m_fluidVelocityGrid.velocityGridU().inBounds(i,j))
            {
                m_fluidVelocityGrid.u(i,j) /= uFluidWeights.at(i,j);
                m_airVelocityGrid.u(i,j) /= uAirWeights.at(i,j);
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                m_fluidVelocityGrid.v(i,j) /= vFluidWeights.at(i,j);
                m_airVelocityGrid.v(i,j) /= vAirWeights.at(i,j);
            }
        }
    }
}

void MultiflipSolver::applyBodyForces()
{
    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(m_fluidVelocityGrid.velocityGridU().inBounds(i,j))
            {
                m_fluidVelocityGrid.u(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().x();
                m_airVelocityGrid.u(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().x();
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                m_fluidVelocityGrid.v(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().y();
                m_airVelocityGrid.v(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().y();
            }
        }
    }
}

void MultiflipSolver::countParticles()
{
    m_fluidParticleCounts.fill(0);
    m_airParticleCounts.fill(0);
    for(MarkerParticle& p : m_markerParticles)
    {
        p.testValue = 0.f;
        int i = std::floor(p.position.x());
        int j = std::floor(p.position.y());
        if(p.material == FluidMaterial::FLUID)
        {
            m_fluidParticleCounts.at(i,j) += 1;
        }
        else
        {
            m_airParticleCounts.at(i,j) += 1;
        }
//        if(m_materialGrid.isSolid(i,j))
//        {
//            std::cout << "Particle in solid at " << i << "," << j << '\n';
//            debug() << "Particle in solid at " << i << "," << j;
//        }
    }
}

void MultiflipSolver::reseedParticles()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            int fluidParticleCount = m_fluidParticleCounts.at(i,j);
            int airParticleCount = m_airParticleCounts.at(i,j);
            if(fluidParticleCount > 20)
            {
                std::cout << "too many particles " << fluidParticleCount << " at " << i << ' ' << j;
            }
            //std::cout << particleCount << " at " << i << " , " << j << std::endl;
            int additionalFluidParticles = SimSettings::particlesPerCell() - fluidParticleCount;
            int additionalAirParticles = SimSettings::particlesPerCell() - airParticleCount;
            if(additionalFluidParticles <= 0 && additionalAirParticles <= 0)
            {
                continue;
            }
            if(m_materialGrid.isSource(i,j))
            {
                for(int p = 0; p < additionalFluidParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    //int emitterId = m_emitterId.at(i,j);
                    float viscosity = 0;
                    float conc = 0;
                    float temp = 0;
                    FluidMaterial mat = FluidMaterial::FLUID;
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc,0.f,mat});
                }
            }
            else if (m_materialGrid.isEmpty(i,j))
            {
                for(int p = 0; p < additionalAirParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_airVelocityGrid.velocityAt(pos);
                    float viscosity = 0;
                    float conc = 0;
                    float temp = 0;
                    FluidMaterial mat = FluidMaterial::EMPTY;
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc,0.f,mat});
                }
            }
//            else if(m_fluidSdf.at(i,j) < SimSettings::particleScale() * SimSettings::dx())
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
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isStrictFluid(i,j) || m_materialGrid.isEmpty(i,j))
            {
                FluidMaterial particleMaterial = m_materialGrid.at(i,j);
                for(int p = 0; p < SimSettings::particlesPerCell(); p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    float viscosity = 0;
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,0.f,0.f,0.f,particleMaterial});
                }
            }
        }
    }
}

void MultiflipSolver::updateMaterials()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
//            if(m_fluidParticleCounts.at(i,j) != 0)
//            {
//                if(m_materialGrid.isEmpty(i,j))
//                {
//                    m_materialGrid.setAt(i,j,FluidMaterial::FLUID);
//                }
//            }
//            else
//            {
//                FluidMaterial m = m_materialGrid.getAt(i,j);
//                if(m == FluidMaterial::FLUID)
//                {
//                    m_materialGrid.setAt(i,j,FluidMaterial::EMPTY);
//                }
//            }
            if(!m_materialGrid.isSolid(i,j))
            {
                if(m_fluidSdf.at(i,j) < 0.f)
                {
                    m_materialGrid.setAt(i,j,FluidMaterial::FLUID);
                }
                else
                {
                    m_materialGrid.setAt(i,j,FluidMaterial::EMPTY);
                }
            }
        }
    }
}

void MultiflipSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double preScale = SimSettings::stepDt() / (SimSettings::dx());

    Grid2d<float> curvatureGrid = simmath::calculateCenteredGridCurvature(m_fluidSdf);

    float minFluidFraction = 0.01;
    float maxFluidFraction = 0.09;

    for (int i = m_sizeI - 1; i >= 0; i--)
    {
        for (int j = m_sizeJ - 1; j >= 0; j--)
        {
            int currentIdx = linearIndex(i,j);
            int im1Idx = linearIndex(i-1,j);
            int jm1Idx = linearIndex(i,j-1);
            Vertex uSamplePos(0.5f + i, j);
            Vertex vSamplePos(i, 0.5f + j);
            double surfaceTensionFactorU = SimSettings::surfaceTensionFactor()
                                            * curvatureGrid.interpolateAt(uSamplePos);
            double surfaceTensionFactorV = SimSettings::surfaceTensionFactor()
                                            * curvatureGrid.interpolateAt(vSamplePos);
            //m_testGrid.at(i,j) = pressures[currentIdx]/100.f;
            //U part
            //m_testGrid.at(i,j) = getWeightedDensityForUSample(i,j);
            if(!m_materialGrid.isSolid(i-1,j) || !m_materialGrid.isSolid(i,j))
            {
                if(m_materialGrid.isSolid(i-1,j) || m_materialGrid.isSolid(i,j))
                {
                    m_fluidVelocityGrid.setU(i,j,0);//Solids are stationary
                    m_airVelocityGrid.setU(i,j,0);
                    m_fluidVelocityGrid.setUValidity(i,j,true);
                    m_airVelocityGrid.setUValidity(i,j,true);
                }
                else
                {
                    float fluidFaceFraction = getFaceFractionUSample(i,j);
                    float invWeightedDensity = 1.f/getWeightedDensityForUSample(i,j);
                    float pressureGrad = pressures[currentIdx]
                                        - (im1Idx != -1? pressures[im1Idx] : 0.f);
                    //m_testGrid.at(i,j) = pressures[currentIdx] / 100000.f;
                    //m_testGrid.at(i,j) = pressureGrad;
                    if(fluidFaceFraction > minFluidFraction)
                    {
                        if(fluidFaceFraction < maxFluidFraction)
                        {//Update both
                            m_fluidVelocityGrid.u(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            m_airVelocityGrid.u(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            if(m_fluidSdf.getAt(i,j) * m_fluidSdf.getAt(i-1,j) < 0.f)
                            {
                                m_fluidVelocityGrid.u(i,j) += preScale * invWeightedDensity * surfaceTensionFactorU;
                                m_airVelocityGrid.u(i,j) += preScale * invWeightedDensity * surfaceTensionFactorU;
                                m_testGrid.at(i,j) = 1.f;
                            }
                            m_fluidVelocityGrid.setUValidity(i,j,true);
                            m_airVelocityGrid.setUValidity(i,j,true);
                        }
                        else
                        {//Only fluid
                            m_fluidVelocityGrid.u(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            if(m_fluidSdf.getAt(i,j) * m_fluidSdf.getAt(i-1,j) < 0.f)
                            {
                                m_fluidVelocityGrid.u(i,j) += preScale * invWeightedDensity * surfaceTensionFactorU;
                            }
                            m_fluidVelocityGrid.setUValidity(i,j,true);
                            m_airVelocityGrid.setUValidity(i,j,false);
                        }
                    }
                    else
                    {//Only air
                        m_airVelocityGrid.u(i,j) -= preScale * invWeightedDensity * pressureGrad;
                        if(m_fluidSdf.getAt(i,j) * m_fluidSdf.getAt(i-1,j) < 0.f)
                        {
                            m_airVelocityGrid.u(i,j) += preScale * invWeightedDensity * surfaceTensionFactorU;
                        }
                        m_fluidVelocityGrid.setUValidity(i,j,false);
                        m_airVelocityGrid.setUValidity(i,j,true);
                    }
                }
            }
            else
            {
                m_fluidVelocityGrid.setUValidity(i,j,false);
                m_airVelocityGrid.setUValidity(i,j,false);
            }

            //V part
            if(!m_materialGrid.isSolid(i,j-1) || !m_materialGrid.isSolid(i,j))
            {
                if(m_materialGrid.isSolid(i,j-1) || m_materialGrid.isSolid(i,j))
                {
                    m_fluidVelocityGrid.setV(i,j,0);//Solids are stationary
                    m_airVelocityGrid.setV(i,j,0);
                    m_fluidVelocityGrid.setVValidity(i,j,true);
                    m_airVelocityGrid.setVValidity(i,j,true);
                }
                else
                {
                    float fluidFaceFraction = getFaceFractionVSample(i,j);
                    float invWeightedDensity = 1.f/getWeightedDensityForVSample(i,j);
                    float pressureGrad = pressures[currentIdx]
                                        - (jm1Idx != -1? pressures[jm1Idx] : 0.f);
                    if(fluidFaceFraction > minFluidFraction)
                    {
                        if(fluidFaceFraction < maxFluidFraction)
                        {//Update both
                            m_fluidVelocityGrid.v(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            m_airVelocityGrid.v(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            if(m_fluidSdf.getAt(i,j) * m_fluidSdf.getAt(i,j-1) < 0.f)
                            {
                                m_fluidVelocityGrid.v(i,j) += preScale * invWeightedDensity * surfaceTensionFactorV;
                                m_airVelocityGrid.v(i,j) += preScale * invWeightedDensity * surfaceTensionFactorV;
                            }
                            m_fluidVelocityGrid.setVValidity(i,j,true);
                            m_airVelocityGrid.setVValidity(i,j,true);
                        }
                        else
                        {//Only fluid
                            m_fluidVelocityGrid.v(i,j) -= preScale * invWeightedDensity * pressureGrad;
                            if(m_fluidSdf.getAt(i,j) * m_fluidSdf.getAt(i,j-1) < 0.f)
                            {
                                m_fluidVelocityGrid.v(i,j) += preScale * invWeightedDensity * surfaceTensionFactorV;
                            }
                            m_fluidVelocityGrid.setVValidity(i,j,true);
                            m_airVelocityGrid.setVValidity(i,j,false);
                        }
                    }
                    else
                    {//Only air
                        m_airVelocityGrid.v(i,j) -= preScale * invWeightedDensity * pressureGrad;
                        if(m_fluidSdf.getAt(i,j) * m_fluidSdf.getAt(i,j-1) < 0.f)
                        {
                            m_airVelocityGrid.v(i,j) += preScale * invWeightedDensity * surfaceTensionFactorV;
                        }
                        m_fluidVelocityGrid.setVValidity(i,j,false);
                        m_airVelocityGrid.setVValidity(i,j,true);
                    }
                }
            }
            else
            {
                m_fluidVelocityGrid.setVValidity(i,j,false);
                m_airVelocityGrid.setVValidity(i,j,false);
            }
        }
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U fluid!\n" << std::flush;
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V fluid!\n" << std::flush;
    }

    if(anyNanInf(m_airVelocityGrid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U air!\n" << std::flush;
    }

    if(anyNanInf(m_airVelocityGrid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V air!\n" << std::flush;
    }
}

void MultiflipSolver::particleUpdate()
{
    Grid2d<float>& prevFluidU = m_savedFluidVelocityGrid.velocityGridU();
    Grid2d<float>& prevFluidV = m_savedFluidVelocityGrid.velocityGridV();
    Grid2d<float>& prevAirU = m_savedAirVelocityGrid.velocityGridU();
    Grid2d<float>& prevAirV = m_savedAirVelocityGrid.velocityGridV();

    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
    {
        MarkerParticle &p = m_markerParticles[i];
        if(p.material == FluidMaterial::FLUID)
        {
            Vertex oldFluidVelocity(simmath::cubicIterpUGrid(p.position.x(),p.position.y(),prevFluidU)
                                    / SimSettings::dx(),
                               simmath::cubicIterpVGrid(p.position.x(),p.position.y(),prevFluidV)
                                    / SimSettings::dx());
            Vertex newFluidVelocity = m_fluidVelocityGrid.velocityAt(p.position) / SimSettings::dx();
            p.velocity = SimSettings::picRatio() * newFluidVelocity +
                    (1.f-SimSettings::picRatio()) * (p.velocity + newFluidVelocity - oldFluidVelocity);
        }
        else
        {
            Vertex oldAirVelocity(simmath::cubicIterpUGrid(p.position.x(),p.position.y(),prevAirU)
                                  / SimSettings::dx(),
                               simmath::cubicIterpVGrid(p.position.x(),p.position.y(),prevAirV)
                                  / SimSettings::dx());
            Vertex newAirVelocity = m_airVelocityGrid.velocityAt(p.position) / SimSettings::dx();
            p.velocity = SimSettings::picRatio() * newAirVelocity +
                    (1.f-SimSettings::picRatio()) * (p.velocity + newAirVelocity - oldAirVelocity);
        }
    }
}

void MultiflipSolver::bumpParticles()
{
    const float escapeRadius = 1.5f * SimSettings::dx();
    const float particleRadius = SimSettings::particleScale() * SimSettings::dx();
    Grid2d<float> curvatureGrid = simmath::calculateCenteredGridCurvature(m_fluidSdf);
    for(MarkerParticle& p : m_markerParticles)
    {
        float sign = p.material == FluidMaterial::FLUID ? -1.f : 1.f;
        float sdfAtParticle = simmath::lerpCenteredGrid(p.position.x(),p.position.y(),m_fluidSdf);
        Vertex sdfGradAtParticle = simmath::gradCenteredGrid(p.position.x(),p.position.y(),m_fluidSdf);
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

            Vertex vectorToClosestSurfacePoint = p.position - m_fluidSdf.closestSurfacePoint(p.position);
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

double MultiflipSolver::divergenceAt(int i, int j)
{
    double weightedUCurrent = getWeightedVelocityUSample(i,j);
    double weightedVCurrent = getWeightedVelocityVSample(i,j);
    double weightedUIp1 = getWeightedVelocityUSample(i+1,j);
    double weightedVJp1 = getWeightedVelocityVSample(i,j+1);

    return (weightedUIp1-weightedUCurrent
            +weightedVJp1-weightedVCurrent);
}

float MultiflipSolver::getFaceFractionUSample(int i, int j)
{
    float sdfCurrent = m_fluidSdf.getAt(i,j);
    float sdfAtIm1 = m_fluidSdf.getAt(i-1,j);

    auto calcD = [this,i,j](float sdfCurrent, float sdfneg)
    {
        float dxSqrd = SimSettings::dx() * SimSettings::dx();
        float deltaSqrd = (sdfCurrent - sdfneg) * (sdfCurrent - sdfneg);
        if(dxSqrd < deltaSqrd)
        {
            return 0.000001f;
        }
        float result = sqrt(dxSqrd - deltaSqrd);
        if(std::isnan(result) || std::isinf(result))
        {
            std::cout << "Nan U d calculation!\n";
        }
        return result;
    };

    float d = calcD(sdfCurrent,sdfAtIm1);
    float result = std::clamp(0.5f - (sdfAtIm1 + sdfCurrent) / (2.f * d),0.f,1.f);
//    if(result < 0.1f)
//    {
//        return 0;
//    }
    if(std::isnan(result) || std::isinf(result))
    {
        std::cout << "Bad U face fraction!\n";
    }
    //m_testGrid.at(i,j) = result;
    return result;
}

float MultiflipSolver::getFaceFractionVSample(int i, int j)
{
    float sdfCurrent = m_fluidSdf.getAt(i,j);
    float sdfAtJm1 = m_fluidSdf.getAt(i,j-1);

    auto calcD = [this,i,j](float sdfCurrent, float sdfneg)
    {
        float dxSqrd = SimSettings::dx() * SimSettings::dx();
        float deltaSqrd = (sdfCurrent - sdfneg) * (sdfCurrent - sdfneg);
        if(dxSqrd < deltaSqrd)
        {
            return 0.000001f;
        }
        float result = sqrt(dxSqrd - deltaSqrd);
        if(std::isnan(result) || std::isinf(result))
        {
            std::cout << "Nan V d calculation!\n";
        }
        return result;
    };

    float d = calcD(sdfCurrent,sdfAtJm1);
    float result = std::clamp(0.5f - (sdfAtJm1 + sdfCurrent) / (2.f * d), 0.f, 1.f);
//    if(result < 0.1f)
//    {
//        return 0;
//    }
    if(std::isnan(result) || std::isinf(result))
    {
        std::cout << "Bad V face fraction!\n";
    }
    //m_testGrid.at(i,j) = result;
    return result;
}

double MultiflipSolver::getWeightedDensityForUSample(int i, int j)
{
    float sdfCurrent = m_fluidSdf.getAt(i,j) / SimSettings::dx();
    float sdfAtIm1 = m_fluidSdf.getAt(i-1,j) / SimSettings::dx();
    Vertex sdfGradAtSample = simmath::gradCenteredGrid(static_cast<float>(i) + 0.5f,
                                                       static_cast<float>(j) + 0.5f,
                                                       m_fluidSdf);
    float normDistCurr = sdfCurrent / sdfGradAtSample.distFromZero();

    Vertex sdfGradAtIm1 = simmath::gradCenteredGrid(static_cast<float>(i) - 0.5f,
                                                       static_cast<float>(j) + 0.5f,
                                                       m_fluidSdf);
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
    fluidFrac = getFaceFractionUSample(i,j);
    fluidFrac = std::clamp(fluidFrac,0.f,1.f);
    float result = SimSettings::fluidDensity() * fluidFrac + (1.f - fluidFrac) * SimSettings::airDensity();
    return result;
}

double MultiflipSolver::getWeightedDensityForVSample(int i, int j)
{
    float sdfCurrent = m_fluidSdf.getAt(i,j) / SimSettings::dx();
    float sdfAtJm1 = m_fluidSdf.getAt(i,j-1) / SimSettings::dx();
    Vertex sdfGradAtSample = simmath::gradCenteredGrid(static_cast<float>(i) + 0.5f,
                                                       static_cast<float>(j) + 0.5f,
                                                       m_fluidSdf);
    float normDistCurr = sdfCurrent / sdfGradAtSample.distFromZero();

    Vertex sdfGradAtJm1 = simmath::gradCenteredGrid(static_cast<float>(i) + 0.5f,
                                                       static_cast<float>(j) - 0.5f,
                                                       m_fluidSdf);
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
    fluidFrac = getFaceFractionVSample(i,j);
    fluidFrac = std::clamp(fluidFrac,0.f,1.f);
    float result = SimSettings::fluidDensity() * fluidFrac + (1.f - fluidFrac) * SimSettings::airDensity();
    return result;
}

double MultiflipSolver::getWeightedVelocityUSample(int i, int j)
{
    double faceFraction = getFaceFractionUSample(i,j);
    double fluidVelocity = m_fluidVelocityGrid.u(i,j);
    double airVelocity = m_airVelocityGrid.u(i,j);

    return fluidVelocity * faceFraction + airVelocity * (1.f - faceFraction);
}

double MultiflipSolver::getWeightedVelocityVSample(int i, int j)
{
    double faceFraction = getFaceFractionVSample(i,j);
    double fluidVelocity = m_fluidVelocityGrid.v(i,j);
    double airVelocity = m_airVelocityGrid.v(i,j);

    return fluidVelocity * faceFraction + airVelocity * (1.f - faceFraction);
}
