#include "flipsolver2d.h"

#include <cmath>
#include <limits>
#include <queue>

#include "linearindexable2d.h"
#include "materialgrid.h"
#include "simsettings.h"
#include "dynamicuppertriangularsparsematrix.h"
#include "logger.h"
#include "mathfuncs.h"

FlipSolver::FlipSolver(int sizeI, int sizeJ, int extrapNeighborRadius, bool vonNeumannNeighbors) :
    LinearIndexable2d(sizeI, sizeJ),
    m_extrapolationNeighborRadius(extrapNeighborRadius),
    m_useVonNeumannNeighborhood(vonNeumannNeighbors),
    m_frameNumber(0),
    m_validVVelocitySampleCount(0),
    m_validUVelocitySampleCount(0),
    m_fluidVelocityGrid(sizeI, sizeJ),
    m_savedFluidVelocityGrid(sizeI, sizeJ),
    m_materialGrid(sizeI,sizeJ, FluidMaterial::SINK),
    m_solidSdf(sizeI,sizeJ),
    m_fluidSdf(sizeI, sizeJ),
    m_knownCenteredParams(sizeI,sizeJ, false, OOBStrategy::OOB_CONST, true),
    m_viscosityGrid(sizeI,sizeJ,0.f,OOBStrategy::OOB_EXTEND),
    m_emitterId(sizeI, sizeJ, -1),
    m_solidId(sizeI, sizeJ,-1),
    m_fluidParticleCounts(sizeI, sizeJ),
    m_divergenceControl(sizeI,sizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f),
    m_testGrid(sizeI, sizeJ)
{
    m_randEngine = std::mt19937(SimSettings::randomSeed());
}

void FlipSolver::setExrapolationRadius(int radius)
{
    m_extrapolationNeighborRadius = radius;
}

void FlipSolver::project()
{
    std::vector<double> rhs(cellCount(),0.0);
    std::vector<double> pressures(cellCount(),0.0);
    calcPressureRhs(rhs);
    //debug() << "Calculated rhs: " << rhs;
    DynamicUpperTriangularSparseMatrix mat = getPressureProjectionMatrix();
    if(!m_pcgSolver.solve(mat,pressures,rhs,SimSettings::pcgIterLimit()))
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
    updateLinearFluidViscosityMapping();
    std::vector<double> rhs(m_validUVelocitySampleCount + m_validVVelocitySampleCount, 0.0);
    std::vector<double> result(m_validUVelocitySampleCount + m_validVVelocitySampleCount, 0.0);
    calcViscosityRhs(rhs);
    DynamicUpperTriangularSparseMatrix mat = getViscosityMatrix();
    //debug() << "mat=" << mat;
    //debug() << "vec=" << rhs;
    //binDump(mat,"test.bin");
    if(!PCGSolver::solve(mat,result,rhs,2000))
    {
        std::cout << "PCG Solver Viscosity: Iteration limit exhaustion!\n";
    }

    if(anyNanInf(result))
    {
        std::cout << "NaN or inf in viscosity vector!\n" << std::flush;
    }

    int vBaseIndex = m_validUVelocitySampleCount;
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            int linearIdxU = linearViscosityVelocitySampleIndexU(i,j);
            int linearIdxV = linearViscosityVelocitySampleIndexV(i,j);
            if(linearIdxU != -1)
            {
                m_fluidVelocityGrid.setU(i,j,result[linearIdxU]);
            }

            if(linearIdxV != -1)
            {
                m_fluidVelocityGrid.setV(i,j,result[vBaseIndex + linearIdxV]);
            }
        }
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U velocity after viscosity!\n" << std::flush;
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridV().data()))
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
        int substepCount = 0;
        p.position = rk3Integrate(p.position,SimSettings::stepDt(), m_fluidVelocityGrid);
        if(m_solidSdf.interpolateAt(p.position.x(),p.position.y()) < 0.f)
        {
            p.position = m_solidSdf.closestSurfacePoint(p.position);
        }
        maxSubsteps = std::max(substepCount,maxSubsteps);
        int pI = simmath::integr(p.position.x());
        int pJ = simmath::integr(p.position.y());
        if(!inBounds(pI,pJ) || m_materialGrid.isSink(pI,pJ))
        {
            m_markerParticles.erase(markerParticles().begin() + i);
        }

    }
    //std::cout << "Advection done in max " << maxSubsteps << " substeps" << std::endl;
}

void FlipSolver::particleUpdate()
{
    Grid2d<float>& prevU = m_savedFluidVelocityGrid.velocityGridU();
    Grid2d<float>& prevV = m_savedFluidVelocityGrid.velocityGridV();
    for(int i = m_markerParticles.size() - 1; i >= 0; i--)
    {
        MarkerParticle &p = m_markerParticles[i];
        Vertex oldVelocity(prevU.interpolateAt(p.position.x(),p.position.y()),
                           prevV.interpolateAt(p.position.x(),p.position.y()));
        Vertex newVelocity = m_fluidVelocityGrid.velocityAt(p.position) ;
//        if(oldVelocity.distFromZero() > (SimSettings::cflNumber() / SimSettings::stepDt()))
//        {
//            p.velocity = newVelocity;
//        }
//        else
//        {
            p.velocity = SimSettings::picRatio() * newVelocity +
                    (1.f-SimSettings::picRatio()) * (p.velocity + newVelocity - oldVelocity);
//            p.temperature = SimSettings::ambientTemp() +
//                    (p.temperature - SimSettings::ambientTemp()) *
//                    std::exp(-SimSettings::tempDecayRate() * SimSettings::stepDt());
//            p.smokeConcentrartion *= std::exp(-SimSettings::concentrartionDecayRate() * SimSettings::stepDt());
//        }
    }
}

void FlipSolver::afterTransfer()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isSource(i,j))
            {
                int emitterId = m_emitterId.at(i,j);
                m_viscosityGrid.at(i,j) = m_sources[emitterId].viscosity();
                if(m_sources[emitterId].velocityTransfer())
                {
                    m_fluidVelocityGrid.setU(i,j,m_sources[emitterId].velocity().x() / SimSettings::dx());
                    m_fluidVelocityGrid.setV(i,j,m_sources[emitterId].velocity().y() / SimSettings::dx());
                    m_fluidVelocityGrid.setUValidity(i,j,true);
                    m_fluidVelocityGrid.setVValidity(i,j,true);
                }
            }
        }
    }
}

void FlipSolver::step()
{
    advect();
    particleToGrid();
    updateSdf();
    updateLinearFluidViscosityMapping();
    updateMaterials();
    afterTransfer();
    extrapolateLevelsetInside(m_fluidSdf);
    m_fluidVelocityGrid.extrapolate(10);

    m_savedFluidVelocityGrid = m_fluidVelocityGrid;
    applyBodyForces();
    project();
    updateVelocityFromSolids();
    applyViscosity();
    project();
    m_fluidVelocityGrid.extrapolate(10);
    particleUpdate();
    countParticles();
    reseedParticles();
}

void FlipSolver::stepFrame()
{
    if(m_frameNumber == 0)
    {
        firstFrameInit();
    }
    m_testGrid.fill(0.f);
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
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
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
            m_solidSdf.at(i,j) = dist;
            if(dist < 0)
            {
                m_materialGrid.at(i,j) = FluidMaterial::SOLID;
                m_solidId.at(i,j) = minIdx;
            }
        }
    }
}

void FlipSolver::updateSources()
{
    float dx = SimSettings::dx();
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            for(int emitterIdx = 0; emitterIdx < m_sources.size(); emitterIdx++)
            {
                Emitter& e = m_sources[emitterIdx];
                if(e.geometry().signedDistance(i*dx,j*dx) <= 0.f)
                {
                    m_materialGrid.at(i,j) = FluidMaterial::SOURCE;
                    m_emitterId.at(i,j) = emitterIdx;
                    m_divergenceControl.at(i,j) = e.divergence();
                }
            }
        }
    }
}

void FlipSolver::updateSinks()
{
    float dx = SimSettings::dx();
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            for(Sink& s : m_sinks)
            {
                if(s.geo().signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) <= 0.f)
                {
                    m_materialGrid.at(i,j) = FluidMaterial::SINK;
                    m_divergenceControl.at(i,j) = s.divergence();
                }
            }
        }
    }
}

void FlipSolver::updateInitialFluid()
{
    float dx = SimSettings::dx();
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            for(Emitter& e : m_initialFluid)
            {
                if(e.geometry().signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx) <= 0.f)
                {
                    m_materialGrid.at(i,j) = FluidMaterial::FLUID;
                    m_viscosityGrid.at(i,j) = e.viscosity();
                }
            }
        }
    }
}

int FlipSolver::gridSizeI()
{
    return m_sizeI;
}

int FlipSolver::gridSizeJ()
{
    return m_sizeJ;
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
    p.velocity = m_fluidVelocityGrid.velocityAt(particle.x(), particle.y());
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

void FlipSolver::reseedParticles()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            int particleCount = m_fluidParticleCounts.at(i,j);
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
            if(m_materialGrid.isSource(i,j))
            {
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    //Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    Vertex velocity = Vertex();
                    int emitterId = m_emitterId.at(i,j);
                    float viscosity = m_sources[emitterId].viscosity();
                    float conc = m_sources[emitterId].concentrartion();
                    float temp = m_sources[emitterId].temperature();
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity,temp,conc});
                }
            }
            else if(m_fluidSdf.at(i,j) < -1.f)
            {
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    float viscosity = m_viscosityGrid.interpolateAt(pos);
                    addMarkerParticle(MarkerParticle{pos,velocity,viscosity});
                }
            }
        }
    }
}

void FlipSolver::seedInitialFluid()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isStrictFluid(i,j))
            {
                for(int p = 0; p < SimSettings::particlesPerCell(); p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    float viscosity = m_viscosityGrid.interpolateAt(pos);
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

const MaterialGrid &FlipSolver::materialGrid() const
{
    return m_materialGrid;
}

const StaggeredVelocityGrid &FlipSolver::fluidVelocityGrid() const
{
    return m_fluidVelocityGrid;
}

const SdfGrid &FlipSolver::fluidSdf() const
{
    return m_fluidSdf;
}

const SdfGrid &FlipSolver::solidSdf() const
{
    return m_solidSdf;
}

const Grid2d<float> &FlipSolver::testGrid() const
{
    return m_testGrid;
}

double FlipSolver::divergenceAt(int i, int j)
{
    return static_cast<double>(m_fluidVelocityGrid.u(i+1,j) - m_fluidVelocityGrid.u(i,j)
                               + m_fluidVelocityGrid.v(i,j+1) - m_fluidVelocityGrid.v(i,j)
                               + m_divergenceControl.at(i,j));
}

std::vector<int> FlipSolver::validSolidNeighborIds(int i, int j)
{
    std::vector<int> output;
    output.reserve(8);
    for (int iOffset = -1; iOffset < 2; iOffset++)
    {
        for (int jOffset = -1; jOffset < 2; jOffset++)
        {
            if(m_solidId.inBounds(i + iOffset, j + jOffset))
            {
                int id = m_solidId.at(i + iOffset, j + jOffset);
                if(id != -1)
                {
                    output.push_back(m_solidId.at(i + iOffset, j + jOffset));
                }
            }
        }
    }

    return output;
}

void FlipSolver::firstFrameInit()
{
    updateInitialFluid();
    seedInitialFluid();
}

void FlipSolver::extrapolateVelocityField(Grid2d<float> &extrapGrid, Grid2d<bool> &flagGrid, int steps)
{
    Grid2d<int> markers(m_sizeI,m_sizeJ,std::numeric_limits<int>().max());
    std::queue<Index2d> wavefront;
    //Extrapolate U
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(flagGrid.at(i,j))
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(Index2d& neighborIndex : getNeighborhood(i,j,m_extrapolationNeighborRadius,m_useVonNeumannNeighborhood))
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
        std::vector<Index2d> neighbors = getNeighborhood(index,m_extrapolationNeighborRadius,m_useVonNeumannNeighborhood);
        double avg = 0;
        int count = 0;
        for(Index2d& neighborIndex : neighbors)
        {
            if(markers.at(neighborIndex) < markers.at(index))
            {
                avg += extrapGrid.at(neighborIndex);
                count++;
            }
            if(markers.at(neighborIndex) == std::numeric_limits<int>().max() && markers.at(index) <= steps)
            {
                markers.at(neighborIndex) = markers.at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        extrapGrid.at(index) = avg / count;
        flagGrid.at(index) = true;

        wavefront.pop();
    }
}

DynamicUpperTriangularSparseMatrix FlipSolver::getPressureProjectionMatrix()
{
    DynamicUpperTriangularSparseMatrix output = DynamicUpperTriangularSparseMatrix(cellCount(),7);

    double scale = SimSettings::stepDt() / (SimSettings::fluidDensity() * SimSettings::dx() * SimSettings::dx());

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    for(int i = 0; i <  m_sizeI; i++)
    {
        for(int j = 0; j <  m_sizeJ; j++)
        {
            if(m_materialGrid.isFluid(i,j))
            {
                //X Neighbors
                if( m_materialGrid.isFluid(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }else if(m_materialGrid.isEmpty(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                if( m_materialGrid.isFluid(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                    output.setAx(i,j,-scale,  indexer);
                } else if(m_materialGrid.isEmpty(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                //Y Neighbors
                if( m_materialGrid.isFluid(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }else if(m_materialGrid.isEmpty(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                if( m_materialGrid.isFluid(i,j+1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                    output.setAy(i,j,-scale,  indexer);
                } else if(m_materialGrid.isEmpty(i,j+1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }
            }
        }
    }

    return output;
}

DynamicUpperTriangularSparseMatrix FlipSolver::getViscosityMatrix()
{
    int validUSamples = m_validUVelocitySampleCount;
    int validVSamples = m_validVVelocitySampleCount;
    DynamicUpperTriangularSparseMatrix output(validUSamples + validVSamples,7);

    float scaleTwoDt = 2*SimSettings::stepDt() / (SimSettings::dx() * SimSettings::dx());
    float scaleTwoDx = SimSettings::stepDt() / (2 * SimSettings::dx() * SimSettings::dx());

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            int currLinearIdxU = linearViscosityVelocitySampleIndexU(i,j);
            int currLinearIdxV = linearViscosityVelocitySampleIndexV(i,j);
            if(currLinearIdxU != -1)
            {
                float fi = static_cast<float>(i);
                float fj = static_cast<float>(j);
                int vBaseIndex = validUSamples;
                //U component
                output.addTo(currLinearIdxU,currLinearIdxU,SimSettings::fluidDensity());

                int uImOneLinearIdx = linearViscosityVelocitySampleIndexU(i-1,j);

                if(uImOneLinearIdx != -1)
                {
                    output.addTo(currLinearIdxU,
                                 uImOneLinearIdx,
                                 -scaleTwoDt * m_viscosityGrid.getAt(i-1,j));

                    output.addTo(currLinearIdxU,
                                 currLinearIdxU,
                                 scaleTwoDt * m_viscosityGrid.getAt(i-1,j));
                }

                int uIpOneLinearIdx = linearViscosityVelocitySampleIndexU(i+1,j);

                if(uIpOneLinearIdx != -1)
                {
                    output.addTo(currLinearIdxU,
                                 uIpOneLinearIdx,
                                 -scaleTwoDt * m_viscosityGrid.getAt(i,j));

                    output.addTo(currLinearIdxU,
                                 currLinearIdxU,
                                 scaleTwoDt * m_viscosityGrid.getAt(i,j));
                }

                int uJmOneLinearIdx = linearViscosityVelocitySampleIndexU(i,j-1);
                int vImOneLinearIdx = linearViscosityVelocitySampleIndexV(i-1,j);

                if(uJmOneLinearIdx != -1
                        && currLinearIdxV != -1
                        && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = m_viscosityGrid.interpolateAt(fi-0.5f,fj-0.5f);
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

                int uJpOneLinearIdx = linearViscosityVelocitySampleIndexU(i,j+1);
                int vJpOneLinearIdx = linearViscosityVelocitySampleIndexV(i,j+1);
                int vImOneJpOneLinearIdx = linearViscosityVelocitySampleIndexV(i-1,j+1);

                if(uJpOneLinearIdx != -1
                        && vJpOneLinearIdx != -1
                        && vImOneJpOneLinearIdx != -1)
                {
                    float lerpedViscosity = m_viscosityGrid.interpolateAt(fi-0.5f,fj+0.5f);
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
                output.addTo(vBaseIndex + currLinearIdxV,vBaseIndex + currLinearIdxV,SimSettings::fluidDensity());

                int vJmOneLinearIdx = linearViscosityVelocitySampleIndexV(i,j-1);

                if(vJmOneLinearIdx != -1)
                {
                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vJmOneLinearIdx,
                                 -scaleTwoDt * m_viscosityGrid.getAt(i,j-1));

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDt * m_viscosityGrid.getAt(i,j-1));
                }

                int vJpOneLinearIdx = linearViscosityVelocitySampleIndexV(i,j+1);

                if(vJpOneLinearIdx != -1)
                {
                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + vJpOneLinearIdx,
                                 -scaleTwoDt * m_viscosityGrid.getAt(i,j));

                    output.addTo(vBaseIndex + currLinearIdxV,
                                 vBaseIndex + currLinearIdxV,
                                 scaleTwoDt * m_viscosityGrid.getAt(i,j));
                }

                int uJmOneLinearIdx = linearViscosityVelocitySampleIndexU(i,j-1);
                int vImOneLinearIdx = linearViscosityVelocitySampleIndexV(i-1,j);

                if(currLinearIdxU != -1
                        && uJmOneLinearIdx != -1
                        && vImOneLinearIdx != -1)
                {
                    float lerpedViscosity = m_viscosityGrid.interpolateAt(fi-0.5f,fj-0.5f);
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

                int uIpOneLinearIdx = linearViscosityVelocitySampleIndexU(i+1,j);
                int uIpOneJmOneLinearIdx = linearViscosityVelocitySampleIndexU(i+1,j-1);
                int vIpOneLinearIdx = linearViscosityVelocitySampleIndexV(i+1,j);

                if(uIpOneLinearIdx != -1
                        && uIpOneJmOneLinearIdx != -1
                        && vIpOneLinearIdx != -1)
                {
                    float lerpedViscosity = m_viscosityGrid.interpolateAt(fi+0.5f,fj-0.5f);
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

int FlipSolver::particleCount()
{
    return m_markerParticles.size();
}

int FlipSolver::cellCount()
{
    return m_sizeI * m_sizeJ;
}

void FlipSolver::calcPressureRhs(std::vector<double> &rhs)
{
    double scale = 1.f/SimSettings::dx();

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if (m_materialGrid.isFluid(i,j))
            {
                rhs[linearIndex(i,j)] = -scale * divergenceAt(i,j);

                if(m_materialGrid.isSolid(i-1,j))
                {
                    rhs[linearIndex(i,j)] -= scale * static_cast<double>(m_fluidVelocityGrid.u(i,j) - 0.f);
                }
                if(m_materialGrid.isSolid(i+1,j))
                {
                    rhs[linearIndex(i,j)] += scale * static_cast<double>(m_fluidVelocityGrid.u(i+1,j) - 0.f);
                }

                if(m_materialGrid.isSolid(i,j-1))
                {
                    rhs[linearIndex(i,j)] -= scale * static_cast<double>(m_fluidVelocityGrid.v(i,j) - 0.f);
                }
                if(m_materialGrid.isSolid(i,j+1))
                {
                    rhs[linearIndex(i,j)] += scale * static_cast<double>(m_fluidVelocityGrid.v(i,j+1) - 0.f);
                }
            }
        }
    }
}

void FlipSolver::calcViscosityRhs(std::vector<double> &rhs)
{
    int vBaseIndex = m_validUVelocitySampleCount;
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            int linearIdxU = linearViscosityVelocitySampleIndexU(i,j);
            int linearIdxV = linearViscosityVelocitySampleIndexV(i,j);
            if(linearIdxU != -1)
            {
                float u = m_fluidVelocityGrid.getU(i,j);
                rhs[linearIdxU] = SimSettings::fluidDensity() * u;
            }

            if(linearIdxV != -1)
            {
                float v = m_fluidVelocityGrid.getV(i,j);
                rhs[vBaseIndex + linearIdxV] = SimSettings::fluidDensity() * v;
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

void FlipSolver::countParticles()
{
    m_fluidParticleCounts.fill(0);
    for(MarkerParticle& p : m_markerParticles)
    {
        int i = std::floor(p.position.x());
        int j = std::floor(p.position.y());
        m_fluidParticleCounts.at(i,j) += 1;
//        if(m_materialGrid.isSolid(i,j))
//        {
//            std::cout << "Particle in solid at " << i << "," << j << '\n';
//            debug() << "Particle in solid at " << i << "," << j;
//        }
    }
}

void FlipSolver::updateMaterials()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_fluidSdf.at(i,j) < 0.f)
            {
                if(m_materialGrid.isEmpty(i,j))
                {
                    m_materialGrid.at(i,j) = FluidMaterial::FLUID;
                }
            }
            else
            {
                if(m_materialGrid.isStrictFluid(i,j))
                {
                    m_materialGrid.at(i,j) = FluidMaterial::EMPTY;
                }
            }
        }
    }
}

void FlipSolver::updateVelocityFromSolids()
{
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            std::vector<int> ids = validSolidNeighborIds(i,j);
            if(!ids.empty())
            {
                float avg = 0;
                for(int& id : ids)
                {
                    Obstacle& obj = m_obstacles[id];
                    avg += obj.friction();
                }
                avg /= ids.size();
                if(m_materialGrid.uSampleAffectedBySolid(i,j))
                {
                    m_fluidVelocityGrid.setU(i,j,m_fluidVelocityGrid.u(i,j) * (1-avg));
                }
                if(m_materialGrid.vSampleAffectedBySolid(i,j))
                {
                    m_fluidVelocityGrid.setV(i,j,m_fluidVelocityGrid.v(i,j) * (1-avg));
                }
            }
        }
    }
}

void FlipSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double scale = SimSettings::stepDt() / (SimSettings::fluidDensity() * SimSettings::dx());

    for (int i = m_sizeI - 1; i >= 0; i--)
    {
        for (int j = m_sizeJ - 1; j >= 0; j--)
        {
            int fluidIndex = linearIndex(i,j);
            int fluidIndexIM1 = linearIndex(i-1,j);
            int fluidIndexJM1 = linearIndex(i,j-1);
            double pCurrent = fluidIndex == -1 ? 0.0 : pressures[fluidIndex];
            //U part
            if(m_materialGrid.isFluid(i-1,j) || m_materialGrid.isFluid(i,j))
            {
                if(m_materialGrid.isSolid(i-1,j) || m_materialGrid.isSolid(i,j))
                {
                    m_fluidVelocityGrid.setU(i,j,0);//Solids are stationary
                }
                else
                {
                    double pIm1 = fluidIndexIM1 == -1 ? 0.0 : pressures[fluidIndexIM1];
                    m_fluidVelocityGrid.u(i,j) -= scale * (pCurrent - pIm1);
                }
            }
            else
            {
                m_fluidVelocityGrid.setUValidity(i,j,false);
            }

            //V part
            if(m_materialGrid.isFluid(i,j-1) || m_materialGrid.isFluid(i,j))
            {
                if(m_materialGrid.isSolid(i,j-1) || m_materialGrid.isSolid(i,j))
                {
                    m_fluidVelocityGrid.setV(i,j,0);//Solids are stationary
                }
                else
                {
                    double pJm1 = fluidIndexJM1 == -1 ? 0.0 : pressures[fluidIndexJM1];
                    m_fluidVelocityGrid.v(i,j) -= scale * (pCurrent - pJm1);
                }
            }
            else
            {
                m_fluidVelocityGrid.setVValidity(i,j,false);
            }
        }
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridU().data()))
    {
        std::cout << "NaN or inf in U vector!\n" << std::flush;
    }

    if(anyNanInf(m_fluidVelocityGrid.velocityGridV().data()))
    {
        std::cout << "NaN or inf in V vector!\n" << std::flush;
    }
}

Vertex FlipSolver::rk3Integrate(Vertex currentPosition, float dt, StaggeredVelocityGrid &grid)
{
    Vertex k1 = grid.velocityAt(currentPosition);
    Vertex k2 = grid.velocityAt(currentPosition + 1.f/2.f * dt * k1);
    Vertex k3 = grid.velocityAt(currentPosition + 3.f/4.f * dt * k2);

    return currentPosition + (2.f/9.f) * dt * k1 + (3.f/9.f) * dt * k2 + (4.f/9.f) * dt * k3;
}

void FlipSolver::particleToGrid()
{
    particleVelocityToGrid();
    centeredParamsToGrid();
}

void FlipSolver::applyBodyForces()
{
    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(m_fluidVelocityGrid.velocityGridU().inBounds(i,j))
            {
                m_fluidVelocityGrid.u(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().x();
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                m_fluidVelocityGrid.v(i,j) +=
                        SimSettings::stepDt() * SimSettings::globalAcceleration().y();
            }
        }
    }
}

void FlipSolver::updateSdf()
{
    float particleRadius = SimSettings::particleScale();
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            float distSqrd = std::numeric_limits<float>::max();
            Vertex centerPoint = Vertex(static_cast<float>(i) + 0.5f,
                                        static_cast<float>(j) + 0.5f);
            for(MarkerParticle& p : m_markerParticles)
            {
                float diffX = p.position.x() - centerPoint.x();
                float diffY = p.position.y() - centerPoint.y();
                float newDistSqrd = diffX*diffX + diffY * diffY;
                if(newDistSqrd < distSqrd)
                {
                    distSqrd = newDistSqrd;
                }
            }

            m_fluidSdf.at(i,j) = std::sqrt(distSqrd) - particleRadius;
        }
    }
}

void FlipSolver::particleVelocityToGrid()
{
    Grid2d<float> uWeights(m_fluidVelocityGrid.velocityGridU().sizeI(),
                           m_fluidVelocityGrid.velocityGridU().sizeJ(),1e-10f);
    Grid2d<float> vWeights(m_fluidVelocityGrid.velocityGridV().sizeI(),
                           m_fluidVelocityGrid.velocityGridV().sizeJ(),1e-10f);
    m_fluidVelocityGrid.velocityGridU().fill(0.f);
    m_fluidVelocityGrid.velocityGridV().fill(0.f);
    m_fluidVelocityGrid.uSampleValidityGrid().fill(false);
    m_fluidVelocityGrid.vSampleValidityGrid().fill(false);

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
                if(uWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        uWeights.at(iIdx,jIdx) += weightU;
                        m_fluidVelocityGrid.u(iIdx,jIdx) += weightU * (p.velocity.x());
                        m_fluidVelocityGrid.setUValidity(iIdx,jIdx,true);
                    }
                }

                if(vWeights.inBounds(iIdx,jIdx))
                {
                    if(std::abs(weightU) > 1e-9f && std::abs(weightV) > 1e-9f)
                    {
                        vWeights.at(iIdx,jIdx) += weightV;
                        m_fluidVelocityGrid.v(iIdx,jIdx) += weightV * (p.velocity.y());
                        m_fluidVelocityGrid.setVValidity(iIdx,jIdx,true);
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
                m_fluidVelocityGrid.u(i,j) /= uWeights.at(i,j);
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                m_fluidVelocityGrid.v(i,j) /= vWeights.at(i,j);
            }
        }
    }
}

void FlipSolver::centeredParamsToGrid()
{
    Grid2d<float> centeredWeights(m_sizeI,m_sizeJ,1e-10f);

    m_viscosityGrid.fill(0.f);
    m_divergenceControl.fill(0.f);
    m_knownCenteredParams.fill(false);

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
                if(inBounds(iIdx,jIdx))
                {
                    float weightCentered = simmath::quadraticBSpline(p.position.x() - (iIdx),
                                                         p.position.y() - (jIdx));
                    if(std::abs(weightCentered) > 1e-6f)
                    {
                        centeredWeights.at(iIdx,jIdx) += weightCentered;
                        m_viscosityGrid.at(iIdx,jIdx) += weightCentered * p.viscosity;
                        m_knownCenteredParams.at(iIdx,jIdx) = true;
                    }
                }
            }
        }
    }

    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(centeredWeights.inBounds(i,j))
            {
                if(m_knownCenteredParams.at(i,j))
                {
                    m_viscosityGrid.at(i,j) /= centeredWeights.at(i,j);
                }
            }
        }
    }
}

void FlipSolver::extrapolateLevelsetInside(SdfGrid &grid)
{
    int sizeI = grid.sizeI();
    int sizeJ = grid.sizeJ();
    Grid2d<int> markers(sizeI,sizeJ,std::numeric_limits<int>().max());
    std::queue<Index2d> wavefront;
    for(int i = 0; i < sizeI; i++)
    {
        for(int j = 0; j < sizeJ; j++)
        {
            if(grid.at(i,j) > 0.f)
            {
                markers.at(i,j) = 0;
            }
        }
    }

    for(int i = 0; i < sizeI; i++)
    {
        for(int j = 0; j < sizeJ; j++)
        {
            if(markers.at(i,j) != 0)
            {
                for(Index2d& neighborIndex : grid.getNeighborhood(i,j,1,false))
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
        std::vector<Index2d> neighbors = grid.getNeighborhood(index,1,false);
        double avg = 0;
        int count = 0;
        for(Index2d& neighborIndex : neighbors)
        {
            if(markers.at(neighborIndex) < markers.at(index))
            {
                avg += grid.at(neighborIndex);
                count++;
            }
            if(markers.at(neighborIndex) == std::numeric_limits<int>().max() && markers.at(index) <= 1e6)
            {
                markers.at(neighborIndex) = markers.at(index) + 1;
                wavefront.push(neighborIndex);
            }
        }
        grid.at(index) = avg / count - 1.f;

        wavefront.pop();
    }
}

void FlipSolver::updateLinearFluidViscosityMapping()
{
    updateValidULinearMapping();
    updateValidVLinearMapping();
}

void FlipSolver::updateValidULinearMapping()
{
    m_uVelocitySamplesMap.clear();
    int linearIdx = 0;
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.uVelocitySampleInside(i,j))
            {
                m_uVelocitySamplesMap.insert(std::make_pair(std::pair(i,j),linearIdx));
                linearIdx++;
            }
        }
    }

    m_validUVelocitySampleCount = linearIdx;
}

void FlipSolver::updateValidVLinearMapping()
{
    m_vVelocitySamplesMap.clear();
    int linearIdx = 0;
    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.vVelocitySampleInside(i,j))
            {
                m_vVelocitySamplesMap.insert(std::make_pair(std::pair(i,j),linearIdx));
                linearIdx++;
            }
        }
    }

    m_validVVelocitySampleCount = linearIdx;
}

int FlipSolver::linearViscosityVelocitySampleIndexU(int i, int j)
{
    std::unordered_map<std::pair<int,int>,int>::iterator iter = m_uVelocitySamplesMap.find(std::pair<int,int>(i,j));
    if(iter != m_uVelocitySamplesMap.end())
    {
        return iter->second;
    }
    else
    {
        return -1;
    }
}

int FlipSolver::linearViscosityVelocitySampleIndexV(int i, int j)
{
    std::unordered_map<std::pair<int,int>,int>::iterator iter = m_vVelocitySamplesMap.find(std::pair<int,int>(i,j));
    if(iter != m_vVelocitySamplesMap.end())
    {
        return iter->second;
    }
    else
    {
        return -1;
    }
}

float FlipSolver::maxParticleVelocity()
{
    float maxVelocitySqr = std::numeric_limits<float>::min();
    for(MarkerParticle& p : m_markerParticles)
    {
        float velocitySqr = p.velocity.x()*p.velocity.x() + p.velocity.y()*p.velocity.y();
        if(velocitySqr > maxVelocitySqr) maxVelocitySqr = velocitySqr;
    }

    return std::sqrt(maxVelocitySqr);
}
