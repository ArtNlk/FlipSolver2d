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
    m_grid.updateLinearFluidViscosityMapping();
    std::vector<double> rhs(m_grid.cellCount(),0.0);
    calcPressureRhs(rhs);
    //debug() << "Calculated rhs: " << rhs;
    DynamicUpperTriangularSparseMatrix mat = getPressureProjectionMatrix();
    std::vector<double> pressures(m_grid.cellCount(),0.0);
    if(!m_pcgSolver.solve(mat,pressures,rhs,200))
    {
        std::cout << "PCG Solver pressure: Iteration limit exhaustion!\n";
    }

    if(anyNanInf(pressures))
    {
        std::cout << "NaN or inf in pressures vector!\n" << std::flush;
    }

    //debug() << "pressures = " << pressures;

    applyPressuresToVelocityField(pressures);
}

void FlipSolver::applyViscosity()
{
    m_grid.updateLinearFluidViscosityMapping();
    int uValidSamples = m_grid.validVelocitySampleCountU();
    int vValidSamples = m_grid.validVelocitySampleCountV();
    std::vector<double> rhs(uValidSamples + vValidSamples, 0);
    std::vector<double> result(uValidSamples + vValidSamples, 0);
    calcViscosityRhs(rhs);
    DynamicUpperTriangularSparseMatrix mat = getViscosityMatrix();
    //debug() << "mat=" << mat;
    //debug() << "vec=" << rhs;
    //binDump(mat,"test.bin");
    if(!m_pcgSolver.solve(mat,result,rhs,200))
    {
        std::cout << "PCG Solver Viscosity: Iteration limit exhaustion!\n";
    }

    if(anyNanInf(result))
    {
        std::cout << "NaN or inf in viscosity vector!\n" << std::flush;
    }

    int vBaseIndex = m_grid.validVelocitySampleCountU();
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            int linearIdxU = m_grid.linearViscosityVelocitySampleIndexU(i,j);
            int linearIdxV = m_grid.linearViscosityVelocitySampleIndexV(i,j);
            if(linearIdxU != -1)
            {
                m_grid.setU(i,j,result[linearIdxU]);
            }

            if(linearIdxV != -1)
            {
                m_grid.setV(i,j,result[vBaseIndex + linearIdxV]);
            }
        }
    }

    if(anyNanInf(m_grid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U velocity after viscosity!\n" << std::flush;
    }

    if(anyNanInf(m_grid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V velocity after viscosity!\n" << std::flush;
    }
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

void FlipSolver::particleUpdate(Grid2d<float> &prevU, Grid2d<float> &prevV)
{
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
            p.velocity = SimSettings::picRatio() * newVelocity +
                    (1.f-SimSettings::picRatio()) * (p.velocity + newVelocity - oldVelocity);
            p.temperature = SimSettings::ambientTemp() +
                    (p.temperature - SimSettings::ambientTemp()) *
                    std::exp(-SimSettings::tempDecayRate() * SimSettings::stepDt());
            p.smokeConcentrartion *= std::exp(-SimSettings::concentrartionDecayRate() * SimSettings::stepDt());
//        }
    }
}

void FlipSolver::step()
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
        particleToGrid();
        break;
    case STAGE_EXTRAPOLATE_VELOCITY:
        extrapolateVelocityField(1);
        break;
    case STAGE_APPLY_GLOBAL_ACCELERATION:
        prevU = m_grid.velocityGridU();
        prevV = m_grid.velocityGridV();
        applyBodyForces();
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

void FlipSolver::updateSolids()
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            float dx = SimSettings::dx();
            float dist = std::numeric_limits<float>::max();
            int minIdx = 0;
            for(int solidIdx = 0; solidIdx < m_obstacles.size(); solidIdx++)
            {
                float sdf = m_obstacles[solidIdx].geometry().signedDistance(
                            (static_cast<float>(i)+0.5)*dx,
                            (static_cast<float>(j)+0.5)*dx) / dx;
                if(sdf < dist)
                {
                    minIdx = solidIdx;
                    dist = sdf;
                }
            }
            m_grid.setSdf(i,j,dist);
            if(dist < 0)
            {
                m_grid.setMaterial(i,j,FluidCellMaterial::SOLID);
                m_grid.setSolidId(i,j,minIdx);
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
            for(int emitterIdx = 0; emitterIdx < m_sources.size(); emitterIdx++)
            {
                Emitter& e = m_sources[emitterIdx];
                if(e.geometry().signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) <= 0.f)
                {
                    m_grid.setMaterial(i,j,FluidCellMaterial::SOURCE);
                    m_grid.setEmitterId(i,j,emitterIdx);
                    m_grid.divergenceControl(i,j) = e.divergence();
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
            for(Sink& s : m_sinks)
            {
                if(s.geo().signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) <= 0.f)
                {
                    m_grid.setMaterial(i,j,FluidCellMaterial::SINK);
                    m_grid.divergenceControl(i,j) = s.divergence();
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
            for(Emitter& e : m_initialFluid)
            {
                if(e.geometry().signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) <= 0.f)
                {
                    m_grid.setMaterial(i,j,FluidCellMaterial::FLUID);
                    m_grid.setViscosity(i,j,e.viscosity());
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

void FlipSolver::addGeometry(Obstacle& geometry)
{
    m_obstacles.push_back(geometry);
}

void FlipSolver::addSource(Emitter &emitter)
{
    m_sources.push_back(emitter);
}

void FlipSolver::addSink(Sink &geometry)
{
    m_sinks.push_back(geometry);
}

void FlipSolver::addInitialFluid(Emitter &emitter)
{
    m_initialFluid.push_back(emitter);
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
                if(additionalParticles <= 0)
                {
                    continue;
                }
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_grid.velocityAt(pos);
                    int emitterId = m_grid.emitterId(i,j);
                    float viscosity = m_sources[emitterId].viscosity();
                    float conc = m_sources[emitterId].concentrartion();
                    float temp = m_sources[emitterId].temperature();
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc});
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
            if(m_grid.isStrictFluid(i,j))
            {
                for(int p = 0; p < SimSettings::particlesPerCell(); p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_grid.velocityAt(pos);
                    float viscosity = m_grid.viscosityAt(pos);
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity});
                }
            }
        }
    }
}

std::vector<Obstacle> &FlipSolver::geometryObjects()
{
    return m_obstacles;
}

std::vector<Emitter> &FlipSolver::sourceObjects()
{
    return m_sources;
}

std::vector<Sink> &FlipSolver::sinkObjects()
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

    markers.fill(std::numeric_limits<int>().max());

    //Extrapolate centered parameters (color, viscosity, etc.)
    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.knownFlagsCenteredParams().at(i,j))
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
                avg += m_grid.viscosityGrid().at(neighborIndex);
                count++;
            }
            if(markers.at(neighborIndex) == std::numeric_limits<int>().max() && markers.at(index) <= steps)
            {
                markers.at(neighborIndex) = markers.at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        m_grid.viscosityGrid().at(index) = avg / count;
        m_grid.knownFlagsCenteredParams().at(index) = true;

        wavefront.pop();
    }
}

DynamicUpperTriangularSparseMatrix FlipSolver::getPressureProjectionMatrix()
{
    DynamicUpperTriangularSparseMatrix output = DynamicUpperTriangularSparseMatrix(m_grid.cellCount(),7);

    double scale = SimSettings::stepDt() / (SimSettings::density() * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i <  m_grid.sizeI(); i++)
    {
        for(int j = 0; j <  m_grid.sizeJ(); j++)
        {
            if( m_grid.isFluid(i,j))
            {
                //X Neighbors
                if( m_grid.isFluid(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }else if( m_grid.isEmpty(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }

                if( m_grid.isFluid(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                    output.setAx(i,j,-scale,  m_grid);
                } else if( m_grid.isEmpty(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }

                //Y Neighbors
                if( m_grid.isFluid(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }else if( m_grid.isEmpty(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }

                if( m_grid.isFluid(i,j+1))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                    output.setAy(i,j,-scale,  m_grid);
                } else if( m_grid.isEmpty(i,j+1))
                {
                    output.addToAdiag(i,j,scale,  m_grid);
                }
            }
        }
    }

    return output;
}

DynamicUpperTriangularSparseMatrix FlipSolver::getViscosityMatrix()
{
    int validUSamples = m_grid.validVelocitySampleCountU();
    int validVSamples = m_grid.validVelocitySampleCountV();
    DynamicUpperTriangularSparseMatrix output(validUSamples + validVSamples,7);

    float scaleTwoDt = 2*SimSettings::stepDt() / (SimSettings::dx() * SimSettings::dx());
    float scaleTwoDx = SimSettings::stepDt() / (2 * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i < m_grid.sizeI(); i++)
    {
        for(int j = 0; j < m_grid.sizeJ(); j++)
        {
            int currLinearIdxU = m_grid.linearViscosityVelocitySampleIndexU(i,j);
            int currLinearIdxV = m_grid.linearViscosityVelocitySampleIndexV(i,j);
            if(currLinearIdxU != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = validUSamples;
                ///swap uip/vjp indexes with current index
                //U component
                output.addTo(currLinearIdxU,currLinearIdxU,SimSettings::density());

                int uImOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexU(i-1,j);

                if(uImOneLinearIdx != -1)
                {
                    output.addTo(currLinearIdxU,
                                 uImOneLinearIdx,
                                 -scaleTwoDt * m_grid.viscosity(i-1,j));

                    output.addTo(currLinearIdxU,
                                 currLinearIdxU,
                                 scaleTwoDt * m_grid.viscosity(i-1,j));
                }

                int uIpOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexU(i+1,j);

                if(uIpOneLinearIdx != -1)
                {
                    output.addTo(currLinearIdxU,
                                 uIpOneLinearIdx,
                                 -scaleTwoDt * m_grid.viscosity(i,j));

                    output.addTo(currLinearIdxU,
                                 currLinearIdxU,
                                 scaleTwoDt * m_grid.viscosity(i,j));
                }

                int uJmOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexU(i,j-1);
                int vImOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexV(i-1,j);

                if(uJmOneLinearIdx != -1
                        && currLinearIdxV != -1
                        && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj-0.5f,m_grid.viscosityGrid());
                    //lerpedViscosity = tempVisc;
                    output.addTo(currLinearIdxU,
                                 uJmOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,
                                 vBaseIndex + vImOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,currLinearIdxU,scaleTwoDx * lerpedViscosity);
                }

                int uJpOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexU(i,j+1);
                int vJpOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexV(i,j+1);
                int vImOneJpOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexV(i-1,j+1);

                if(uJpOneLinearIdx != -1
                        && vJpOneLinearIdx != -1
                        && vImOneJpOneLinearIdx != -1)
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj+0.5f,m_grid.viscosityGrid());
                    //lerpedViscosity = tempVisc;
                    output.addTo(currLinearIdxU,
                                 uJpOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,
                                 vBaseIndex + vJpOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,
                                 vBaseIndex + vImOneJpOneLinearIdx,
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(currLinearIdxU,currLinearIdxU,scaleTwoDx * lerpedViscosity);
                }
            }

            //V component
            if(currLinearIdxV != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = validUSamples;
                output.addTo(vBaseIndex + currLinearIdxV,vBaseIndex + currLinearIdxV,SimSettings::density());

                int vJmOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexV(i,j-1);

                if(vJmOneLinearIdx != -1)
                {
                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vJmOneLinearIdx,
                                 -scaleTwoDt * m_grid.viscosity(i,j-1));

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDt * m_grid.viscosity(i,j-1));
                }

                int vJpOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexV(i,j+1);

                if(vJpOneLinearIdx != -1)
                {
                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vJpOneLinearIdx,
                                 -scaleTwoDt * m_grid.viscosity(i,j));

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDt * m_grid.viscosity(i,j));
                }

                int uJmOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexU(i,j-1);
                int vImOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexV(i-1,j);

                if(currLinearIdxU != -1
                        && uJmOneLinearIdx != -1
                        && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi-0.5f,fj-0.5f,m_grid.viscosityGrid());
                    //lerpedViscosity = tempVisc;

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 currLinearIdxU,
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 uJmOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vImOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDx * lerpedViscosity);
                }

                int uIpOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexU(i+1,j);
                int uIpOneJmOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexU(i+1,j-1);
                int vIpOneLinearIdx = m_grid.linearViscosityVelocitySampleIndexV(i+1,j);

                if(uIpOneLinearIdx != -1
                        && uIpOneJmOneLinearIdx != -1
                        && vIpOneLinearIdx != -1)
                {
                    float lerpedViscosity = math::lerpCenteredGrid(fi+0.5f,fj-0.5f,m_grid.viscosityGrid());
                    //lerpedViscosity = tempVisc;

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 uIpOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 uIpOneJmOneLinearIdx,
                                 scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vIpOneLinearIdx,
                                 -scaleTwoDx * lerpedViscosity);

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDx * lerpedViscosity);
                }
            }
        }
    }


    return output;
}

void FlipSolver::calcPressureRhs(std::vector<double> &rhs)
{
    double scale = 1/SimSettings::dx();

    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if (m_grid.isFluid(i,j))
            {
                rhs[m_grid.linearIndex(i,j)] = -scale * static_cast<double>(m_grid.getU(i+1,j) - m_grid.getU(i,j)
                                                              +m_grid.getV(i,j+1) - m_grid.getV(i,j));

                if(m_grid.isSolid(i-1,j))
                {
                    rhs[m_grid.linearIndex(i,j)] -= scale * static_cast<double>(m_grid.getU(i,j) - 0);
                }
                if(m_grid.isSolid(i+1,j))
                {
                    rhs[m_grid.linearIndex(i,j)] += scale * static_cast<double>(m_grid.getU(i+1,j) - 0);
                }

                if(m_grid.isSolid(i,j-1))
                {
                    rhs[m_grid.linearIndex(i,j)] -= scale * static_cast<double>(m_grid.getV(i,j) - 0);
                }
                if(m_grid.isSolid(i,j+1))
                {
                    rhs[m_grid.linearIndex(i,j)] += scale * static_cast<double>(m_grid.getV(i,j+1) - 0);
                }
            }
        }
    }
}

void FlipSolver::calcViscosityRhs(std::vector<double> &rhs)
{
    int vBaseIndex = m_grid.validVelocitySampleCountU();
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            int linearIdxU = m_grid.linearViscosityVelocitySampleIndexU(i,j);
            int linearIdxV = m_grid.linearViscosityVelocitySampleIndexV(i,j);
            if(linearIdxU != -1)
            {
                float u = m_grid.getU(i,j);
                rhs[linearIdxU] = SimSettings::density() * u;
            }

            if(linearIdxV != -1)
            {
                float v = m_grid.getV(i,j);
                rhs[vBaseIndex + linearIdxV] = SimSettings::density() * v;
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

void FlipSolver::updateVelocityFromSolids()
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            std::vector<int> ids = m_grid.nearValidSolidIds(i,j);
            if(!ids.empty())
            {
                float avg = 0;
                for(int& id : ids)
                {
                    Obstacle& obj = m_obstacles[id];
                    avg += obj.friction();
                }
                avg /= ids.size();
                if(m_grid.uSampleAffectedBySolid(i,j))
                {
                    m_grid.setU(i,j,m_grid.getU(i,j) * (1-avg));
                }
                if(m_grid.vSampleAffectedBySolid(i,j))
                {
                    m_grid.setV(i,j,m_grid.getV(i,j) * (1-avg));
                }
            }
        }
    }
}

void FlipSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double scale = SimSettings::stepDt() / (SimSettings::density() * SimSettings::dx());

    for (int i = m_grid.sizeI() - 1; i >= 0; i--)
    {
        for (int j = m_grid.sizeJ() - 1; j >= 0; j--)
        {
            int fluidIndex = m_grid.linearIndex(i,j);
            int fluidIndexIM1 = m_grid.linearIndex(i-1,j);
            int fluidIndexJM1 = m_grid.linearIndex(i,j-1);
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

    if(anyNanInf(m_grid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U vector!\n" << std::flush;
    }

    if(anyNanInf(m_grid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V vector!\n" << std::flush;
    }
}

Vertex FlipSolver::rk3Integrate(Vertex currentPosition, float dt)
{
    Vertex k1 = m_grid.velocityAt(currentPosition);
    Vertex k2 = m_grid.velocityAt(currentPosition + 1.f/2.f * dt * k1);
    Vertex k3 = m_grid.velocityAt(currentPosition + 3.f/4.f * dt * k2);

    return currentPosition + (2.f/9.f) * dt * k1 + (3.f/9.f) * dt * k2 + (4.f/9.f) * dt * k3;
}

void FlipSolver::particleToGrid()
{
    ASSERT(m_grid.velocityGridU().sizeI() == m_grid.velocityGridU().sizeI() && m_grid.velocityGridU().sizeJ() == m_grid.velocityGridU().sizeJ());
    ASSERT(m_grid.velocityGridV().sizeI() == m_grid.velocityGridV().sizeI() && m_grid.velocityGridV().sizeJ() == m_grid.velocityGridV().sizeJ());

    Grid2d<float> uWeights(m_grid.velocityGridU().sizeI(),m_grid.velocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vWeights(m_grid.velocityGridV().sizeI(),m_grid.velocityGridV().sizeJ(),1e-10f);
    Grid2d<float> centeredWeights(m_grid.sizeI(),m_grid.sizeJ(),1e-10f);

    m_grid.velocityGridU().fill(0.f);
    m_grid.velocityGridV().fill(0.f);
    m_grid.viscosityGrid().fill(0.f);
    m_grid.knownFlagsGridU().fill(false);
    m_grid.knownFlagsGridV().fill(false);
    m_grid.knownFlagsCenteredParams().fill(false);
    m_grid.smokeConcentrationGrid().fill(0.f);
    m_grid.temperatureGrid().fill(0);

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
                float weightCentered = math::quadraticBSpline(p.position.x() - (iIdx),
                                                     p.position.y() - (jIdx));
                if(uWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        uWeights.at(iIdx,jIdx) += weightU;
                        m_grid.velocityGridU().at(iIdx,jIdx) += weightU * (p.velocity.x() * SimSettings::dx());
                        m_grid.knownFlagsGridU().at(iIdx,jIdx) = true;
                    }
                }

                if(vWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        vWeights.at(iIdx,jIdx) += weightV;
                        m_grid.velocityGridV().at(iIdx,jIdx) += weightV * (p.velocity.y() * SimSettings::dx());
                        m_grid.knownFlagsGridV().at(iIdx,jIdx) = true;
                    }
                }

                if(centeredWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightCentered) > 1e-9f)
                    {
                        centeredWeights.at(iIdx,jIdx) += weightCentered;
                        m_grid.viscosityGrid().at(iIdx,jIdx) += weightCentered * p.viscosity;
                        m_grid.temperatureGrid().at(iIdx,jIdx) += weightCentered * p.temperature;
                        m_grid.smokeConcentrationGrid().at(iIdx,jIdx) += weightCentered * p.smokeConcentrartion;
                        m_grid.knownFlagsCenteredParams().at(iIdx,jIdx) = true;
                    }
                }
            }
        }
    }

    for (int i = 0; i < m_grid.sizeI() + 1; i++)
    {
        for (int j = 0; j < m_grid.sizeJ() + 1; j++)
        {
            if(m_grid.velocityGridU().inBounds(i,j))
            {
                m_grid.velocityGridU().at(i,j) /= uWeights.at(i,j);
            }
            if(m_grid.velocityGridV().inBounds(i,j))
            {
                m_grid.velocityGridV().at(i,j) /= vWeights.at(i,j);
            }
            if(m_grid.knownFlagsCenteredParams().inBounds(i,j))
            {
                if(m_grid.knownFlagsCenteredParams().at(i,j))
                {
                    m_grid.viscosityGrid().at(i,j) /= centeredWeights.at(i,j);
                    m_grid.smokeConcentrationGrid().at(i,j) /= centeredWeights.at(i,j);
                    m_grid.temperatureGrid().at(i,j) /= centeredWeights.at(i,j);
                }
                else
                {
                    m_grid.temperatureGrid().at(i,j) = SimSettings::ambientTemp();
                }
            }
        }
    }
}

void FlipSolver::applyBodyForces()
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
