#include "flipsolver2d.h"

#include <limits>
#include <queue>

#include "simsettings.h"
#include "dynamicuppertriangularsparsematrix.h"
#include "logger.h"
#include "mathfuncs.h"

FlipSolver::FlipSolver(int extrapRadius, bool vonNeumannNeighbors) :
    m_extrapolationRadius(extrapRadius),
    m_useVonNeumannNeighborhood(vonNeumannNeighbors),
    m_frameNumber(0),
    m_grid(SimSettings::gridSizeI(), SimSettings::gridSizeJ()),
    m_stepStage(SimulationStepStage::STAGE_ADVECT)
{
    m_randEngine = std::mt19937(SimSettings::randomSeed());
}

void FlipSolver::init()
{
    m_grid = MACFluidGrid(SimSettings::gridSizeI(), SimSettings::gridSizeJ());
}

void FlipSolver::project()
{
    m_grid.updateLinearFluidCellMapping();
    std::vector<double> rhs(m_grid.fluidCellCount(),0.0);
    calcRhs(rhs);
    debug() << "Calculated rhs: " << rhs;
    DynamicUpperTriangularSparseMatrix dmat = DynamicUpperTriangularSparseMatrix::forPressureProjection(m_grid);
    UpperTriangularMatrix mat = UpperTriangularMatrix(dmat);
    std::vector<double> pressures(m_grid.fluidCellCount(),0.0);
    if(!m_pcgSolver.solve(mat,m_grid,pressures,rhs,200))
    {
        std::cout << "PCG Solver: Iteration limit exhaustion!\n";
    }

    debug() << "pressures = " << pressures;

    double scale = SimSettings::stepDt() / (SimSettings::density() * SimSettings::dx());

    for (int i = m_grid.sizeI() - 1; i >= 0; i--)
    {
        for (int j = m_grid.sizeJ() - 1; j >= 0; j--)
        {
            int fluidIndex = m_grid.linearFluidIndex(i,j);
            int fluidIndexIM1 = m_grid.linearFluidIndex(i-1,j);
            int fluidIndexJM1 = m_grid.linearFluidIndex(i,j-1);
            //U part
            if(m_grid.isFluid(i-1,j) || m_grid.isFluid(i,j))
            {
                if(m_grid.isSolid(i-1,j) || m_grid.isSolid(i,j))
                {
                    m_grid.setU(i,j,0);//Solids are stationary
                }
                else
                {
                    m_grid.velocityGridU().at(i,j) -= scale * (fluidIndex == -1 ? 0.0 : pressures[fluidIndex] - (fluidIndexIM1 == -1 ? 0.0 : pressures[fluidIndexIM1]));
                }
            }
            else
            {
                m_grid.knownFlagsGridU().at(i,j) = false;
            }

            //V part
            if(m_grid.isFluid(i,j-1) || m_grid.isFluid(i,j))
            {
                if(m_grid.isSolid(i,j-1) || m_grid.isSolid(i,j))
                {
                    m_grid.velocityGridV().at(i,j) = 0;//Solids are stationary
                }
                else
                {
                    m_grid.velocityGridV().at(i,j) -= scale * (fluidIndex == -1 ? 0.0 : pressures[fluidIndex] - (fluidIndexJM1 == -1 ? 0.0 : pressures[fluidIndexJM1]));
                }
            }
            else
            {
                m_grid.knownFlagsGridV().at(i,j) = false;
            }
        }
    }
}

void FlipSolver::applyViscosity()
{
    std::vector<double> currentVelocities(m_grid.fluidCellCount() * 2);
    std::vector<double> newVelocities(m_grid.fluidCellCount() * 2,0);
    m_grid.getFlattenedFluidVelocities(currentVelocities);
    DynamicUpperTriangularSparseMatrix dmat = DynamicUpperTriangularSparseMatrix::forViscosity(m_grid);
    UpperTriangularMatrix mat = UpperTriangularMatrix(dmat);
    if(!m_pcgSolver.solve(mat,m_grid,newVelocities,currentVelocities,200))
    {
        std::cout << "PCG Solver Viscosity: Iteration limit exhaustion!\n";
    }

    m_grid.unflattenFluidVelocities(newVelocities);
}

void FlipSolver::advect()
{
    int maxSubsteps = -1;
    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
    {
        MarkerParticle &p = m_markerParticles[i];
//        float substepTime = 0.f;
//        bool finished = false;
        int substepCount = 0;
//        while(!finished)
//        {
//            Vertex velocity = m_grid.velocityAt(p.position.x(),p.position.y());
//            float maxSubstepSize = 1.f/(velocity.distFromZero() + 1e-15f);
//            if(substepTime + maxSubstepSize >= SimSettings::stepDt())
//            {
//                maxSubstepSize = SimSettings::stepDt() - substepTime;
//                finished = true;
//            }
//            else if(substepTime + 2.f*maxSubstepSize >= SimSettings::stepDt())
//            {
//                maxSubstepSize = 0.5f*(SimSettings::stepDt() - substepTime);
//            }
//            substepTime += maxSubstepSize;
//            substepCount++;
//        }
        p.position = rk3Integrate(p.position,SimSettings::stepDt());
        if(m_grid.sdfAt(p.position.x(),p.position.y()) < 0.f)
        {
            p.position = m_grid.closestSurfacePoint(p.position);
        }
        maxSubsteps = std::max(substepCount,maxSubsteps);
        int pI = math::integr(p.position.x());
        int pJ = math::integr(p.position.y());
        if(!m_grid.inBounds(pI,pJ) || m_grid.isSink(pI,pJ))
        {
            m_markerParticles.erase(markerParticles().begin() + i);
        }

    }
    std::cout << "Advection done in max " << maxSubsteps << " substeps" << std::endl;
}

void FlipSolver::step()
{
    Grid2d<int> particleCounts(m_grid.sizeI(), m_grid.sizeJ());
    m_grid.updateLinearFluidCellMapping();
    countParticles(particleCounts);
    reseedParticles(particleCounts);
    updateMaterialsFromParticles(particleCounts);
    particleVelocityToGrid(m_grid.velocityGridU(),m_grid.velocityGridV());
    extrapolateVelocityField(1);
    Grid2d<float> prevU = m_grid.velocityGridU();
    Grid2d<float> prevV = m_grid.velocityGridV();
    applyGlobalAcceleration();
    project();
    applyViscosity();
    extrapolateVelocityField(1);
    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
    {
        MarkerParticle &p = m_markerParticles[i];
        Vertex oldVelocity(math::lerpUGrid(p.position.x(),p.position.y(),prevU) / SimSettings::dx(),math::lerpVGrid(p.position.x(),p.position.y(),prevV) / SimSettings::dx());
        Vertex newVelocity = m_grid.velocityAt(p.position) / SimSettings::dx();
//        if(oldVelocity.distFromZero() > (SimSettings::cflNumber() / SimSettings::stepDt()))
//        {
//            p.velocity = newVelocity;
//        }
//        else
//        {
            p.velocity = SimSettings::picRatio() * newVelocity + (1.f-SimSettings::picRatio()) * (p.velocity + newVelocity - oldVelocity);
//        }
    }
    advect();
}

void FlipSolver::stagedStep()
{
    static Grid2d<int> particleCounts(m_grid.sizeI(), m_grid.sizeJ());
    static Grid2d<float> prevU(m_grid.velocityGridU().sizeI(), m_grid.velocityGridU().sizeJ());
    static Grid2d<float> prevV(m_grid.velocityGridV().sizeI(), m_grid.velocityGridV().sizeJ());
    m_stepStage++;
    switch(m_stepStage)
    {
    case STAGE_RESEED:
        if(m_frameNumber == 0)
        {
            updateInitialFluid();
            seedInitialFluid();
        }
        particleCounts.fill(0);
        countParticles(particleCounts);
        reseedParticles(particleCounts);
        break;
    case STAGE_UPDATE_MATERIALS:
        updateMaterialsFromParticles(particleCounts);
        break;
    case STAGE_TRANSFER_PARTICLE_VELOCITY:
        particleVelocityToGrid(m_grid.velocityGridU(),m_grid.velocityGridV());
        break;
    case STAGE_EXTRAPOLATE_VELOCITY:
        extrapolateVelocityField(1);
        break;
    case STAGE_APPLY_GLOBAL_ACCELERATION:
        prevU = m_grid.velocityGridU();
        prevV = m_grid.velocityGridV();
        applyGlobalAcceleration();
        break;
    case STAGE_PROJECT:
        project();
        break;
    case STAGE_EXTRAPOLATE_AFTER_PROJECTION:
        extrapolateVelocityField(1);
        break;
    case STAGE_UPDATE_PARTICLE_VELOCITIES:
        for(int i = m_markerParticles.size() - 1; i >= 0; i--)
        {
            MarkerParticle &p = m_markerParticles[i];
            Vertex oldVelocity(math::lerpUGrid(p.position.x(),p.position.y(),prevU) / SimSettings::dx(),math::lerpVGrid(p.position.x(),p.position.y(),prevV) / SimSettings::dx());
            Vertex newVelocity = m_grid.velocityAt(p.position) / SimSettings::dx();
            p.velocity = SimSettings::picRatio() * newVelocity + (1.f-SimSettings::picRatio()) * (p.velocity + newVelocity - oldVelocity);
        }
        break;
    case STAGE_ADVECT:
        advect();
        m_frameNumber++;
        std::cout << "Simulation step finished" << std::endl;
        break;
    case STAGE_ITER_END:
        break;
    }
}

void FlipSolver::stepFrame()
{
    if(m_frameNumber == 0)
    {
        updateInitialFluid();
        seedInitialFluid();
    }
    float substepTime = 0.f;
    bool finished = false;
    int substepCount = 0;
//    for(int i = 0; i < SimSettings::maxSubsteps(); i++)
//    {
//        step();
//    }
    while(!finished)
    {
        float vel = maxParticleVelocity();
        float maxSubstepSize = SimSettings::cflNumber()/(vel + 1e-15f);
        if(substepTime + maxSubstepSize >= SimSettings::frameDt() ||
                substepCount == (SimSettings::maxSubsteps() - 1))
        {
            maxSubstepSize = SimSettings::frameDt() - substepTime;
            finished = true;
        }
        else if(substepTime + 2.f*maxSubstepSize >= SimSettings::frameDt())
        {
            maxSubstepSize = 0.5f*(SimSettings::frameDt() - substepTime);
        }
        SimSettings::stepDt() = maxSubstepSize;
        std::cout << "Substep " << substepCount << " substep dt: " << SimSettings::stepDt() << " vel " << vel << std::endl;
        step();
        substepTime += maxSubstepSize;
        substepCount++;
        if(substepCount > 50) break;
    }
    std::cout << "Frame done in " << substepCount << " substeps" << std::endl;
    m_frameNumber++;
}

void FlipSolver::updateSdf()
{
    float dx = SimSettings::dx();
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            float dist = std::numeric_limits<float>().max();
            for(Geometry2d& geo : m_geometry)
            {
                dist = std::min(geo.signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) / dx,dist);
            }
            m_grid.setSdf(i,j,dist);
        }
    }
}

void FlipSolver::updateSolids()
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.sdf(i,j) < 0)
            {
                m_grid.setMaterial(i,j,FluidCellMaterial::SOLID);
            }
        }
    }
}

void FlipSolver::updateSources()
{
    float dx = SimSettings::dx();
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            for(Geometry2d& geo : m_sources)
            {
                if(geo.signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) <= 0.f)
                {
                    m_grid.setMaterial(i,j,FluidCellMaterial::SOURCE);
                }
            }
        }
    }
}

void FlipSolver::updateSinks()
{
    float dx = SimSettings::dx();
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            for(Geometry2d& geo : m_sinks)
            {
                if(geo.signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) <= 0.f)
                {
                    m_grid.setMaterial(i,j,FluidCellMaterial::SINK);
                }
            }
        }
    }
}

void FlipSolver::updateInitialFluid()
{
    float dx = SimSettings::dx();
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            for(Geometry2d& geo : m_initialFluid)
            {
                if(geo.signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) <= 0.f)
                {
                    m_grid.setMaterial(i,j,FluidCellMaterial::FLUID);
                }
            }
        }
    }
}

int FlipSolver::gridSizeI()
{
    return m_grid.sizeI();
}

int FlipSolver::gridSizeJ()
{
    return m_grid.sizeJ();
}

void FlipSolver::addGeometry(Geometry2d& geometry)
{
    m_geometry.push_back(geometry);
}

void FlipSolver::addSource(Geometry2d &geometry)
{
    m_sources.push_back(geometry);
}

void FlipSolver::addSink(Geometry2d &geometry)
{
    m_sinks.push_back(geometry);
}

void FlipSolver::addInitialFluid(Geometry2d &geometry)
{
    m_initialFluid.push_back(geometry);
}

void FlipSolver::addMarkerParticle(Vertex particle)
{
    MarkerParticle p;
    p.position = particle;
    p.velocity = m_grid.velocityAt(particle.x(), particle.y());
    m_markerParticles.push_back(p);
}

void FlipSolver::addMarkerParticle(MarkerParticle particle)
{
    m_markerParticles.push_back(particle);
}

int FlipSolver::frameNumber()
{
    return m_frameNumber;
}

void FlipSolver::reseedParticles(Grid2d<int> &particleCounts)
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.isSource(i,j))
            {
                int particleCount = particleCounts.at(i,j);
                if(particleCount > 20)
                {
                    std::cout << "too many particles " << particleCount << " at " << i << ' ' << j;
                }
                //std::cout << particleCount << " at " << i << " , " << j << std::endl;
                int additionalParticles = SimSettings::particlesPerCell() - particleCount;
                if(additionalParticles <= 0) continue;
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_grid.velocityAt(pos);
                    addMarkerParticle(MarkerParticle{pos,velocity});
                }
            }
        }
    }
}

void FlipSolver::seedInitialFluid()
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.isFluid(i,j))
            {
                for(int p = 0; p < SimSettings::particlesPerCell(); p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_grid.velocityAt(pos);
                    addMarkerParticle(MarkerParticle{pos,velocity});
                }
            }
        }
    }
}

std::vector<Geometry2d> &FlipSolver::geometryObjects()
{
    return m_geometry;
}

std::vector<Geometry2d> &FlipSolver::sourceObjects()
{
    return m_sources;
}

std::vector<Geometry2d> &FlipSolver::sinkObjects()
{
    return m_sinks;
}

std::vector<MarkerParticle> &FlipSolver::markerParticles()
{
    return m_markerParticles;
}

SimulationStepStage FlipSolver::stepStage()
{
    return m_stepStage;
}

void FlipSolver::extrapolateVelocityField(int steps)
{
    Grid2d<int> markers(m_grid.sizeI(),m_grid.sizeJ(),std::numeric_limits<int>().max());
    std::queue<Index2d> wavefront;
    //Extrapolate U
    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.knownFlagsGridU().at(i,j))
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(Index2d& neighborIndex : m_grid.getNeighborhood(i,j,m_extrapolationRadius,m_useVonNeumannNeighborhood))
                {
                    if(markers.at(neighborIndex) == 0)
                    {
                        markers.at(i,j) = 1;
                        wavefront.push(Index2d(i,j));
                        break;
                    }
                }
            }
        }
    }

    while(!wavefront.empty())
    {
        Index2d index = wavefront.front();
        std::vector<Index2d> neighbors = m_grid.getNeighborhood(index,m_extrapolationRadius,m_useVonNeumannNeighborhood);
        double avg = 0;
        int count = 0;
        for(Index2d& neighborIndex : neighbors)
        {
            if(markers.at(neighborIndex) < markers.at(index))
            {
                avg += m_grid.velocityGridU().at(neighborIndex);
                count++;
            }
            if(markers.at(neighborIndex) == std::numeric_limits<int>().max() && markers.at(index) <= steps)
            {
                markers.at(neighborIndex) = markers.at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        m_grid.velocityGridU().at(index) = avg / count;
        m_grid.knownFlagsGridU().at(index) = true;

        wavefront.pop();
    }

    markers.fill(std::numeric_limits<int>().max());

    //Extrapolate V
    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.knownFlagsGridV().at(i,j))
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(Index2d& neighborIndex : m_grid.getNeighborhood(i,j,m_extrapolationRadius,m_useVonNeumannNeighborhood))
                {
                    if(markers.at(neighborIndex) == 0)
                    {
                        markers.at(i,j) = 1;
                        wavefront.push(Index2d(i,j));
                        break;
                    }
                }
            }
        }
    }

    while(!wavefront.empty())
    {
        Index2d index = wavefront.front();
        std::vector<Index2d> neighbors = m_grid.getNeighborhood(index,m_extrapolationRadius,m_useVonNeumannNeighborhood);
        double avg = 0;
        int count = 0;
        for(Index2d& neighborIndex : neighbors)
        {
            if(markers.at(neighborIndex) < markers.at(index))
            {
                avg += m_grid.velocityGridV().at(neighborIndex);
                count++;
            }
            if(markers.at(neighborIndex) == std::numeric_limits<int>().max() && markers.at(index) <= steps)
            {
                markers.at(neighborIndex) = markers.at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        m_grid.velocityGridV().at(index) = avg / count;
        m_grid.knownFlagsGridV().at(index) = true;

        wavefront.pop();
    }
}

void FlipSolver::calcRhs(std::vector<double> &rhs)
{
    double scale = 1/SimSettings::dx();

    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if (m_grid.isFluid(i,j))
            {
                rhs[m_grid.linearFluidIndex(i,j)] = -scale * static_cast<double>(m_grid.getU(i+1,j) - m_grid.getU(i,j)
                                                              +m_grid.getV(i,j+1) - m_grid.getV(i,j));

                if(m_grid.isSolid(i-1,j))
                {
                    rhs[m_grid.linearFluidIndex(i,j)] -= scale * static_cast<double>(m_grid.getU(i,j) - 0);
                }
                if(m_grid.isSolid(i+1,j))
                {
                    rhs[m_grid.linearFluidIndex(i,j)] += scale * static_cast<double>(m_grid.getU(i+1,j) - 0);
                }

                if(m_grid.isSolid(i,j-1))
                {
                    rhs[m_grid.linearFluidIndex(i,j)] -= scale * static_cast<double>(m_grid.getV(i,j) - 0);
                }
                if(m_grid.isSolid(i,j+1))
                {
                    rhs[m_grid.linearFluidIndex(i,j)] += scale * static_cast<double>(m_grid.getV(i,j+1) - 0);
                }
            }
        }
    }
}

Vertex FlipSolver::jitteredPosInCell(int i, int j)
{
    static std::uniform_real_distribution<float> dist(0.f,1.f);
    float x = (static_cast<float>(i) + dist(m_randEngine));
    float y = (static_cast<float>(j) + dist(m_randEngine));
    return Vertex(x,y);
}

void FlipSolver::countParticles(Grid2d<int> &output)
{
    ASSERT(output.sizeI() == m_grid.sizeI());
    ASSERT(output.sizeJ() == m_grid.sizeJ());
    output.fill(0);
    for(MarkerParticle& p : m_markerParticles)
    {
        int i = std::floor(p.position.x());
        int j = std::floor(p.position.y());
        output.at(i,j) += 1;
//        if(m_grid.isSolid(i,j))
//        {
//            std::cout << "Particle in solid at " << i << "," << j << '\n';
//            debug() << "Particle in solid at " << i << "," << j;
//        }
    }
}

void FlipSolver::updateMaterialsFromParticles(Grid2d<int> &particleCount)
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(particleCount.at(i,j) != 0)
            {
                if(m_grid.isEmpty(i,j))
                {
                    m_grid.setMaterial(i,j,FluidCellMaterial::FLUID);
                }
            }
            else
            {
                FluidCellMaterial m = m_grid.getMaterial(i,j);
                if(m == FluidCellMaterial::FLUID)
                {
                    m_grid.setMaterial(i,j,FluidCellMaterial::EMPTY);
                }
            }
        }
    }
}

Vertex FlipSolver::rk3Integrate(Vertex currentPosition, float dt)
{
    Vertex k1 = m_grid.velocityAt(currentPosition);
    Vertex k2 = m_grid.velocityAt(currentPosition + 1.f/2.f * dt * k1);
    Vertex k3 = m_grid.velocityAt(currentPosition + 3.f/4.f * dt * k2);

    return currentPosition + (2.f/9.f) * dt * k1 + (3.f/9.f) * dt * k2 + (4.f/9.f) * dt * k3;
}

void FlipSolver::particleVelocityToGrid(Grid2d<float> &gridU, Grid2d<float> &gridV)
{
    ASSERT(gridU.sizeI() == m_grid.velocityGridU().sizeI() && gridU.sizeJ() == m_grid.velocityGridU().sizeJ());
    ASSERT(gridV.sizeI() == m_grid.velocityGridV().sizeI() && gridV.sizeJ() == m_grid.velocityGridV().sizeJ());

    Grid2d<float> uWeights(gridU.sizeI(),gridU.sizeJ(),1e-10f);
    Grid2d<float> vWeights(gridV.sizeI(),gridV.sizeJ(),1e-10f);

    gridU.fill(0.f);
    gridV.fill(0.f);
    m_grid.knownFlagsGridU().fill(false);
    m_grid.knownFlagsGridV().fill(false);

    for(MarkerParticle &p : m_markerParticles)
    {
        int i = math::integr(p.position.x());
        int j = math::integr(p.position.y());
        //Run over all cells that this particle might affect
        for (int iOffset = -3; iOffset < 3; iOffset++)
        {
            for (int jOffset = -3; jOffset < 3; jOffset++)
            {
                int iIdx = i+iOffset;
                int jIdx = j+jOffset;
                float weightU = math::quadraticBSpline(p.position.x() - (iIdx),
                                                     p.position.y() - (static_cast<float>(jIdx) + 0.5f));

                float weightV = math::quadraticBSpline(p.position.x() - (static_cast<float>(iIdx) + 0.5f),
                                                     p.position.y() - (jIdx));
                if(uWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        uWeights.at(iIdx,jIdx) += weightU;
                        gridU.at(iIdx,jIdx) += weightU * (p.velocity.x() * SimSettings::dx());
                        m_grid.knownFlagsGridU().at(iIdx,jIdx) = true;
                    }
                }

                if(vWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        vWeights.at(iIdx,jIdx) += weightV;
                        gridV.at(iIdx,jIdx) += weightV * (p.velocity.y() * SimSettings::dx());
                        m_grid.knownFlagsGridV().at(iIdx,jIdx) = true;
                    }
                }
            }
        }
    }

    for (int i = 0; i < m_grid.sizeI() + 1; i++)
    {
        for (int j = 0; j < m_grid.sizeJ() + 1; j++)
        {
            if(gridU.inBounds(i,j)) gridU.at(i,j) /= uWeights.at(i,j);
            if(gridV.inBounds(i,j)) gridV.at(i,j) /= vWeights.at(i,j);
        }
    }
}

void FlipSolver::applyGlobalAcceleration()
{
    for (int i = 0; i < m_grid.sizeI() + 1; i++)
    {
        for (int j = 0; j < m_grid.sizeJ() + 1; j++)
        {
            if(m_grid.velocityGridU().inBounds(i,j)) m_grid.velocityGridU().at(i,j) += SimSettings::stepDt() * SimSettings::globalAcceleration().x();
            if(m_grid.velocityGridV().inBounds(i,j)) m_grid.velocityGridV().at(i,j) += SimSettings::stepDt() * SimSettings::globalAcceleration().y();
        }
    }
}

float FlipSolver::maxParticleVelocity()
{
    float maxVelocitySqr = std::numeric_limits<float>().min();
    for(MarkerParticle& p : m_markerParticles)
    {
        float velocitySqr = p.velocity.x()*p.velocity.x() + p.velocity.y()*p.velocity.y();
        if(velocitySqr > maxVelocitySqr) maxVelocitySqr = velocitySqr;
    }

    return std::sqrt(maxVelocitySqr);
}
