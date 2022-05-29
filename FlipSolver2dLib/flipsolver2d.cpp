#include "flipsolver2d.h"

#include <limits>
#include <queue>

#include "simsettings.h"
#include "dynamicuppertriangularsparsematrix.h"
#include "logger.h"
#include "functions.h"

FlipSolver::FlipSolver(int sizeX, int sizeY, double fluidDensity, double timestepSize,  double sideLength, int extrapRadius, bool vonNeumannNeighbors) :
    m_grid(sizeX, sizeY),
    m_extrapolationRadius(extrapRadius),
    m_useVonNeumannNeighborhood(vonNeumannNeighbors)
{
    SimSettings::density() = fluidDensity;
    SimSettings::dt() = timestepSize;
    SimSettings::dx() = sideLength;
    SimSettings::randomSeed() = 0;
    SimSettings::particlesPerCell() = 8;
    m_randEngine = std::mt19937(SimSettings::randomSeed());
}

void FlipSolver::project()
{
    m_grid.updateLinearToFluidMapping();
    std::vector<double> rhs(m_grid.fluidCellCount(),0.0);
    calcRhs(rhs);
    debug() << "Calculated rhs: " << rhs;
    DynamicUpperTriangularSparseMatrix dmat = DynamicUpperTriangularSparseMatrix(m_grid);
    UpperTriangularMatrix mat = UpperTriangularMatrix(dmat);
    std::vector<double> pressures(m_grid.fluidCellCount(),0.0);
    if(!m_pcgSolver.solve(mat,m_grid,pressures,rhs,200))
    {
        std::cout << "PCG Solver: Iteration limit exhaustion!\n";
    }

    debug() << "pressures = " << pressures;

    double scale = SimSettings::dt() / (SimSettings::density() * SimSettings::dx());

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

void FlipSolver::advect()
{
    Grid2d<int> particleCounts(m_grid.sizeI(), m_grid.sizeJ());
    countParticles(particleCounts);
    reseedParticles(particleCounts);
    updateMaterialsFromParticles(particleCounts);
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
                dist = std::min(geo.signedDistance((static_cast<float>(i)+0.5)*dx,(static_cast<float>(j)+0.5)*dx),dist);
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

void FlipSolver::reseedParticles(Grid2d<int> &particleCounts)
{
    for (int i = 0; i < m_grid.sizeI(); i++)
    {
        for (int j = 0; j < m_grid.sizeJ(); j++)
        {
            if(m_grid.isSource(i,j))
            {
                int particleCount = particleCounts.at(i,j);
                //std::cout << particleCount << " at " << i << " , " << j << std::endl;
                int additionalParticles = SimSettings::particlesPerCell() - particleCount;
                if(additionalParticles <= 0) continue;
                for(int p = 0; p < additionalParticles; p++)
                {
                    addMarkerParticle(jitteredPosInCell(i,j));
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

void FlipSolver::extrapolateVelocityField()
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
            if(markers.at(neighborIndex) == std::numeric_limits<int>().max())
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
            if(markers.at(neighborIndex) == std::numeric_limits<int>().max())
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
                rhs[m_grid.linearFluidIndex(i,j)] = -scale * (m_grid.getU(i+1,j) - m_grid.getU(i,j)
                                                              +m_grid.getV(i,j+1) - m_grid.getV(i,j));

                if(m_grid.isSolid(i-1,j))
                {
                    rhs[m_grid.linearFluidIndex(i,j)] -= scale * (m_grid.getU(i,j) - 0);
                }
                if(m_grid.isSolid(i+1,j))
                {
                    rhs[m_grid.linearFluidIndex(i,j)] += scale * (m_grid.getU(i+1,j) - 0);
                }

                if(m_grid.isSolid(i,j-1))
                {
                    rhs[m_grid.linearFluidIndex(i,j)] -= scale * (m_grid.getV(i,j) - 0);
                }
                if(m_grid.isSolid(i,j+1))
                {
                    rhs[m_grid.linearFluidIndex(i,j)] += scale * (m_grid.getV(i,j+1) - 0);
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
        if(m_grid.isSolid(i,j))
        {
            std::cout << "Particle in solid at " << i << "," << j << '\n';
            debug() << "Particle in solid at " << i << "," << j;
        }
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
