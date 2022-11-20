#include "multiflipsolver.h"
#include "flipsolver2d.h"
#include "fluidcell.h"
#include "geometry2d.h"
#include "grid2d.h"
#include "mathfuncs.h"
#include "simsettings.h"
#include <cmath>
#include <functional>
#include <limits>

MultiflipSolver::MultiflipSolver(int extrapRadius, bool vonNeumannNeighbors) :
    FlipSolver(extrapRadius, vonNeumannNeighbors)
{

}

void MultiflipSolver::step()
{
    updateSdf();
    m_grid.updateLinearFluidViscosityMapping();
    countParticles();
    reseedParticles();
    combineSdf();
    bumpParticles();
    updateMaterialsFromParticles();
    particleToGrid();
    extrapolateVelocityField(1);
    Grid2d<float> prevU = m_grid.velocityGridU();
    Grid2d<float> prevV = m_grid.velocityGridV();
    applyBodyForces();
    project();
    updateVelocityFromSolids();
    //applyViscosity();
    project();
    extrapolateVelocityField(1);
    particleUpdate(prevU, prevV);
    advect();
}

void MultiflipSolver::updateSdf()
{
    float particleRadius = SimSettings::particleScale();
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

            m_grid.fluidSdfGrid().at(i,j) = std::sqrt(fluidDistSqrd * SimSettings::dx()) - particleRadius;
            m_grid.airSdfGrid().at(i,j) = std::sqrt(airDistSqrd * SimSettings::dx()) - particleRadius;
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
            m_grid.fluidSdfGrid().at(i,j) = (m_grid.fluidSdfAt(i,j)-m_grid.airSdfGrid().at(i,j))/2.f;
        }
    }
}

void MultiflipSolver::extrapolateSdf(Grid2d<float> &sdfGrid)
{
    Grid2d<float> normalDeriv(sdfGrid.sizeI(), sdfGrid.sizeJ(), 1e5f);
    Grid2d<bool> extrapFlags(sdfGrid.sizeI(), sdfGrid.sizeJ(), true);

    std::function<float (Grid2d<float> &, Vertex &, void*)> updateFunc = simmath::normalDerivLinearExapolationUpdate;

    for (int i = 0; i < sdfGrid.sizeI(); i++)
    {
        for (int j = 0; j < sdfGrid.sizeJ(); j++)
        {
            if(sdfGrid.at(i,j) > 0)
            {
                Vertex grad = simmath::gradCenteredGrid(i,j,sdfGrid);
                Vertex normal = grad.normalized();
                normalDeriv.at(i,j) = normal.dot(grad);
                extrapFlags.at(i,j) = false;
            }
        }
    }
    simmath::fastSweep(normalDeriv,extrapFlags,updateFunc, nullptr);

    //sdfGrid = normalDeriv;

    updateFunc = simmath::sdfLinearExapolationUpdate;

    simmath::fastSweep(sdfGrid,extrapFlags,updateFunc, &normalDeriv);

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

void MultiflipSolver::countParticles()
{
    m_grid.fluidParticleCountGrid().fill(0);
    m_grid.airParticleCountGrid().fill(0);
    for(MarkerParticle& p : m_markerParticles)
    {
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
            int particleCount = m_grid.fluidParticleCountGrid().at(i,j);
            if(particleCount > 20)
            {
                std::cout << "too many particles " << particleCount << " at " << i << ' ' << j;
            }
            //std::cout << particleCount << " at " << i << " , " << j << std::endl;
            int additionalParticles = SimSettings::particlesPerCell() - particleCount;
            if(additionalParticles <= 0)
            {
                continue;
            }
            if(m_grid.isSource(i,j))
            {
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_grid.fluidVelocityAt(pos);
                    int emitterId = m_grid.emitterId(i,j);
                    float viscosity = m_sources[emitterId].viscosity();
                    float conc = m_sources[emitterId].concentrartion();
                    float temp = m_sources[emitterId].temperature();
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc});
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

void MultiflipSolver::updateMaterialsFromParticles()
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
        }
    }
}

void MultiflipSolver::bumpParticles()
{
    const float escapeRadius = 1.5f;
    Grid2d<float> curvatureGrid = simmath::calculateCenteredGridCurvature(m_grid.fluidSdfGrid());
    for(MarkerParticle& p : m_markerParticles)
    {
        float sign = p.material == FluidMaterial::FLUID ? -1.f : 1.f;
        float sdfAtParticle = simmath::lerpCenteredGrid(p.position.x(),p.position.y(),m_grid.fluidSdfGrid());
        Vertex sdfGradAtParticle = simmath::gradCenteredGrid(p.position.x(),p.position.y(),m_grid.fluidSdfGrid());
        float d = sdfAtParticle / std::abs(sdfGradAtParticle.distFromZero());
        if(std::abs(sdfAtParticle) < escapeRadius)
        {
            float curvature = simmath::lerpCenteredGrid(p.position, curvatureGrid);
            const float curvatureMin = 1.f/(SimSettings::dx()*4);
            const float curvatureMax = 1.f/(SimSettings::dx()*2);
            float negativeSignCurvature = -sign*curvature;
            float targetDist = 0.f;
            if(negativeSignCurvature < curvatureMin)
            {
                targetDist = SimSettings::particleScale();
            }

            if(negativeSignCurvature <= curvatureMax)
            {
                if(sign*d < 0)
                {
                    Vertex vectorToNearestPoint = m_grid.closesFluidSurfacePoint(p.position) - p.position;
                    float distToNearestPoint = vectorToNearestPoint.distFromZero();
                    if(distToNearestPoint > 0)
                    {
                        Vertex newPosition = p.position - (targetDist + distToNearestPoint)*
                                (vectorToNearestPoint/distToNearestPoint);
                        if(std::isnan(newPosition.x()) || std::isnan(newPosition.y()) || std::isnan(newPosition.z()))
                        {
                            std::cout << "Nan in particle bumping";
                        }
                        p.position = newPosition;
                    }
                }
                else if(sign*d < targetDist)
                {
                    Vertex vectorToNearestPoint = m_grid.closesFluidSurfacePoint(p.position) - p.position;
                    float distToNearestPoint = vectorToNearestPoint.distFromZero();
                    if(distToNearestPoint > 0)
                    {
                        Vertex newPosition = p.position + (targetDist - distToNearestPoint)*
                                (vectorToNearestPoint/distToNearestPoint);
                        if(std::isnan(newPosition.x()) || std::isnan(newPosition.y()) || std::isnan(newPosition.z()))
                        {
                            std::cout << "Nan in particle bumping";
                        }
                        p.position = newPosition;
                    }
                }
            }
        }
    }
}
