#include "flipsmokesolver.h"

#include "flipsolver2d.h"
#include "mathfuncs.h"

#include <cmath>

FlipSmokeSolver::FlipSmokeSolver(const SmokeSolverParameters* p):
    FlipSolver(p),
    m_temperature(p->gridSizeI, p->gridSizeJ, p->ambientTemperature, OOBStrategy::OOB_CONST, p->ambientTemperature),
    m_smokeConcentration(p->gridSizeI, p->gridSizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f),
    m_advectedVelocity(p->gridSizeI, p->gridSizeJ),
    m_ambientTemperature(p->ambientTemperature),
    m_temperatureDecayRate(p->temperatureDecayRate),
    m_concentrationDecayRate(p->concentrationDecayRate),
    m_buoyancyFactor(p->buoyancyFactor),
    m_sootFactor(p->sootFactor)
{
    m_projectTolerance = 1e-6;
    m_viscosityEnabled = false;
}

void FlipSmokeSolver::applyBodyForces()
{
    const float sootWeight = m_sootFactor;
    const float buoyancyInfluence = m_buoyancyFactor/m_ambientTemperature;
    const float factor = m_stepDt / m_dx;
    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            if(m_fluidVelocityGrid.velocityGridU().inBounds(i,j))
            {
                float tempAt = m_temperature.lerpolateAt(i,static_cast<float>(j)+0.5f);
                float tempDiff = tempAt - m_ambientTemperature;
                float concentration = m_smokeConcentration.lerpolateAt(i,static_cast<float>(j)+0.5f);
                float accelerationU = (sootWeight * concentration - buoyancyInfluence * tempDiff) *
                                        m_globalAcceleration.x() * factor;
                m_fluidVelocityGrid.u(i,j) += accelerationU;
            }
            if(m_fluidVelocityGrid.velocityGridV().inBounds(i,j))
            {
                float tempAt = m_temperature.lerpolateAt(static_cast<float>(i)+0.5f,j);
                float tempDiff = tempAt - m_ambientTemperature;
                float concentration = m_smokeConcentration.lerpolateAt(static_cast<float>(i)+0.5f,j);
                float accelerationV = (sootWeight * concentration - buoyancyInfluence * tempDiff) *
                                        m_globalAcceleration.y() * factor;
                m_fluidVelocityGrid.v(i,j) += accelerationV;
            }
        }
    }
}

void FlipSmokeSolver::centeredParamsToGrid()
{
    Grid2d<float> centeredWeights(m_sizeI,m_sizeJ,1e-10f);

    m_smokeConcentration.fill(0.f);
    m_temperature.fill(0.f);
    m_divergenceControl.fill(0.f);
    m_knownCenteredParams.fill(false);

    std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
                                              (m_temperatureIndex);

    std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
                                               (m_concentrationIndex);

    for(size_t idx = 0; idx < m_fluidVelocityGrid.linearSize(); idx++)
    {
        Index2d i2d = m_fluidVelocityGrid.index2d(idx);
        std::array<int,9> affectingBins = m_markerParticles.binsForGridCell(i2d);
        for(int binIdx : affectingBins)
        {
            if(binIdx >= 0)
            {
                for(size_t particleIdx : m_markerParticles.binForBinIdx(binIdx))
                {
                    Vertex position = m_markerParticles.particlePosition(particleIdx);
                    float weightCentered = simmath::quadraticBSpline(position.x() - (i2d.i),
                                                                     position.y() - (i2d.j));
                    if(std::abs(weightCentered) > 1e-6f)
                    {
                        centeredWeights.at(i2d.i,i2d.j) += weightCentered;
                        m_temperature.at(i2d.i,i2d.j) +=
                            weightCentered * particleTemperatures.at(particleIdx);

                        m_smokeConcentration.at(i2d.i,i2d.j) +=
                            weightCentered * particleConcentrations.at(particleIdx);

                        m_knownCenteredParams.at(i2d.i,i2d.j) = true;
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
                    m_temperature.at(i,j) /= centeredWeights.at(i,j);
                    m_smokeConcentration.at(i,j) /= centeredWeights.at(i,j);
//                    m_temperature.at(i,j) = m_ambientTemperature +
//                                             (m_temperature.at(i,j) - m_ambientTemperature) *
//                                                 std::exp(-m_temperatureDecayRate * m_stepDt);
                    //m_smokeConcentration.at(i,j) *= std::exp(-m_concentrationDecayRate * m_stepDt);
                }
                else
                {
                    m_temperature.at(i,j) = m_ambientTemperature;
                }
                m_testGrid.at(i,j) = m_temperature.at(i,j) / 1000.f;
            }
        }
    }
}

void FlipSmokeSolver::calcPressureRhs(std::vector<double> &rhs)
{
    double scale = 1/m_dx;

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if (!m_materialGrid.isSolid(i,j))
            {
                rhs[linearIndex(i,j)] = -scale * divergenceAt(i,j);

                if(m_materialGrid.isSolid(i-1,j))
                {
                    rhs[linearIndex(i,j)] -= scale * static_cast<double>(m_fluidVelocityGrid.u(i,j) - 0);
                }
                if(m_materialGrid.isSolid(i+1,j))
                {
                    rhs[linearIndex(i,j)] += scale * static_cast<double>(m_fluidVelocityGrid.u(i+1,j) - 0);
                }

                if(m_materialGrid.isSolid(i,j-1))
                {
                    rhs[linearIndex(i,j)] -= scale * static_cast<double>(m_fluidVelocityGrid.v(i,j) - 0);
                }
                if(m_materialGrid.isSolid(i,j+1))
                {
                    rhs[linearIndex(i,j)] += scale * static_cast<double>(m_fluidVelocityGrid.v(i,j+1) - 0);
                }
            }
        }
    }
}

void FlipSmokeSolver::particleUpdate()
{
    FlipSolver::particleUpdate();
    std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
                                               (m_temperatureIndex);

    std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
                                                 (m_concentrationIndex);

    for(int i = 0; i < m_markerParticles.particleCount(); i++)
    {
        Vertex &position = m_markerParticles.particlePosition(i);
        //particleTemperatures.at(i) = m_temperature.lerpolateAt(position);
        particleTemperatures.at(i) = m_ambientTemperature +
                (particleTemperatures.at(i) - m_ambientTemperature) *
                std::exp(-m_temperatureDecayRate * m_stepDt);
        //particleConcentrations.at(i) = m_smokeConcentration.lerpolateAt(position);
        particleConcentrations.at(i) *= std::exp(-m_concentrationDecayRate * m_stepDt);
        //p.testValue = p.smokeConcentrartion;
    }
}

void FlipSmokeSolver::afterTransfer()
{
    FlipSolver::afterTransfer();
    m_divergenceControl.fill(0.f);
    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            if(m_materialGrid.isSource(i,j))
            {
                int emitterId = m_emitterId.at(i,j);
                m_smokeConcentration.at(i,j) = m_sources[emitterId].concentrartion();
                m_temperature.at(i,j) = m_sources[emitterId].temperature();
            }
        }
    }
}

void FlipSmokeSolver::reseedParticles()
{
    for (int pIndex = 0; pIndex < m_markerParticles.particleCount(); pIndex++)
    {
        Vertex pos = m_markerParticles.particlePosition(pIndex);
        int i = pos.x();
        int j = pos.y();

        if(m_fluidParticleCounts.at(i,j) > 2*m_particlesPerCell)
        {
            m_markerParticles.markForDeath(pIndex);
            m_fluidParticleCounts.at(i,j) -= 1;
            //pIndex--;
            continue;
        }
    }

    std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
                                               (m_temperatureIndex);

    std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
                                                 (m_concentrationIndex);

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
            int additionalParticles = (m_particlesPerCell / 2) - particleCount;
            if(additionalParticles <= 0)
            {
                continue;
            }
            if(m_materialGrid.isSource(i,j))
            {
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
//                    int emitterId = m_emitterId.at(i,j);
//                    //float viscosity = m_sources[emitterId].viscosity();
//                    float conc = m_sources[emitterId].concentrartion();
//                    float temp = m_sources[emitterId].temperature();
                    float conc = m_smokeConcentration.interpolateAt(pos);
                    float temp = m_temperature.interpolateAt(pos);
                    size_t pIdx = m_markerParticles.addMarkerParticle(pos,velocity);
                    particleTemperatures[pIdx] = temp;
                    particleConcentrations[pIdx] = conc;
                    //testValues[pIdx] = temp / 1000.f;
                }
            }
        }
    }
}

void FlipSmokeSolver::applyPressuresToVelocityField(std::vector<double> &pressures)
{
    double scale = m_stepDt / (m_fluidDensity * m_dx);

    for (int i = m_sizeI - 1; i >= 0; i--)
    {
        for (int j = m_sizeJ - 1; j >= 0; j--)
        {
//            if(!m_materialGrid.isFluid(i,j))
//            {
//                continue;
//            }
            int fluidIndex = linearIndex(i,j);
            int fluidIndexIM1 = linearIndex(i-1,j);
            int fluidIndexJM1 = linearIndex(i,j-1);
            double pCurrent = fluidIndex == -1 ? 0.0 : pressures[fluidIndex];
            //U part
            if(!m_materialGrid.isSolid(i-1,j) || !m_materialGrid.isSolid(i,j))
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
            if(!m_materialGrid.isSolid(i,j-1) || !m_materialGrid.isSolid(i,j))
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

void FlipSmokeSolver::step()
{
    advect();
    //densityCorrection();
    m_markerParticles.pruneParticles();
    m_markerParticles.rebinParticles();
    particleToGrid();

    m_savedFluidVelocityGrid = m_fluidVelocityGrid;

    updateSdf();
    //extrapolateLevelsetOutside(m_fluidSdf);
    //updateLinearFluidViscosityMapping();
    updateMaterials();
    afterTransfer();
    m_fluidVelocityGrid.extrapolate(5);
    extrapolateLevelsetOutside(m_fluidSdf);
    combineAdvectedGrids();

    applyBodyForces();

    project();

    updateVelocityFromSolids();

    //m_fluidVelocityGrid.extrapolate(5);
    particleUpdate();
    countParticles();
    reseedParticles();
}

void FlipSmokeSolver::advect()
{
    FlipSolver::advect();

    Vertex offsetU(0.f,0.5f);
    Vertex offsetV(0.5f,0.f);

    Grid2d<float> advectedU(m_sizeI + 1, m_sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetU);
    Grid2d<float> advectedV(m_sizeI, m_sizeJ + 1, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetV);

    std::vector<Range> rangesU = ThreadPool::i()->splitRange(advectedU.data().size());
    std::vector<Range> rangesV = ThreadPool::i()->splitRange(advectedV.data().size());

    for(Range r : rangesU)
    {
        ThreadPool::i()->enqueue(&FlipSmokeSolver::eulerAdvectionThread,this,r,
                                 offsetU,std::ref(m_fluidVelocityGrid.velocityGridU()), std::ref(advectedU));
    }
    for(Range r : rangesV)
    {
        ThreadPool::i()->enqueue(&FlipSmokeSolver::eulerAdvectionThread,this,r,
                                 offsetV,std::ref(m_fluidVelocityGrid.velocityGridV()), std::ref(advectedV));
    }
    ThreadPool::i()->wait();

    m_advectedVelocity.velocityGridU() = advectedU;
    m_advectedVelocity.velocityGridV() = advectedV;
}

LinearSolver::MatElementProvider FlipSmokeSolver::getPressureMatrixElementProvider()
{
    return std::bind(&FlipSmokeSolver::getMatFreeElementForLinIdx,this,std::placeholders::_1);
}

LinearSolver::SparseMatRowElements FlipSmokeSolver::getMatFreeElementForLinIdx(unsigned int i)
{
    double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);
    std::array<int,4> neighbors = immidiateNeighbors(static_cast<int>(i));
    std::vector<FluidMaterial>& materials = m_materialGrid.data();

    LinearSolver::SparseMatRowElements output;
    output.fill(std::pair<int, double>(0,0.0));
    output[4].first = i;

    if(!solidTest(materials[i]))
    {
        for(unsigned int i = 0; i < neighbors.size(); i++)
        {
            if(neighbors[i] < 0 || neighbors[i] >= materials.size())
            {
                continue;
            }
            output[i].first = neighbors[i];
            output[i].second = -scale * fluidTest(materials[neighbors[i]]);
            output[4].second += scale * !solidTest(materials[neighbors[i]]);
        }
    }

    return output;
}

DynamicUpperTriangularSparseMatrix FlipSmokeSolver::getPressureProjectionMatrix()
{
    DynamicUpperTriangularSparseMatrix output = DynamicUpperTriangularSparseMatrix(cellCount(),7);

    double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    for(int i = 0; i <  m_sizeI; i++)
    {
        for(int j = 0; j <  m_sizeJ; j++)
        {
            if(!m_materialGrid.isSolid(i,j))
            {
                //X Neighbors
                if(m_materialGrid.isFluid(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }else if(m_materialGrid.isEmpty(i-1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                if(m_materialGrid.isFluid(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                    output.setAx(i,j,-scale,  indexer);
                } else if(m_materialGrid.isEmpty(i+1,j))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                //Y Neighbors
                if(m_materialGrid.isFluid(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }else if(m_materialGrid.isEmpty(i,j-1))
                {
                    output.addToAdiag(i,j,scale,  indexer);
                }

                if(m_materialGrid.isFluid(i,j+1))
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

void FlipSmokeSolver::eulerAdvectionThread(Range range, Vertex offset, const Grid2d<float> &inputGrid, Grid2d<float> &outputGrid)
{
    std::vector<float>& dataOut = outputGrid.data();
    for(int idx = range.start; idx < range.end; idx++)
    {
        Index2d i2d = outputGrid.index2d(idx);
        Vertex currentPos = Vertex(i2d.i, i2d.j) + offset;
        Vertex prevPos = inverseRk4Integrate(currentPos,m_fluidVelocityGrid);
//        if(m_solidSdf.interpolateAt(prevPos) < 0.f)
//        {
//            prevPos = m_solidSdf.closestSurfacePoint(prevPos);
//        }
        dataOut[idx] = inputGrid.interpolateAt(prevPos);
    }
}

Vertex FlipSmokeSolver::inverseRk4Integrate(Vertex newPosition, StaggeredVelocityGrid &grid)
{
    float factor = -m_stepDt;
    Vertex k1 = factor*grid.velocityAt(newPosition);
    Vertex k2 = factor*grid.velocityAt(newPosition + 0.5f*k1);
    Vertex k3 = factor*grid.velocityAt(newPosition + 0.5f*k2);
    Vertex k4 = factor*grid.velocityAt(newPosition + k3);

    return newPosition + (1.0f/6.0f)*(k1 + 2.f*k2 + 2.f*k3 + k4);
}

void FlipSmokeSolver::combineAdvectedGrids()
{
    auto nbcombine = [this](float vAdvected, float vParticle, float sdf, Vertex pos)
    {
        //if(sdf < 1.f || m_solidSdf.interpolateAt(pos) < 2.f)
        if(true)
        {
            return vParticle;
        }
        else
        {
            return vAdvected;
        }
    };

    for (int i = 0; i < m_sizeI + 1; i++)
    {
        for (int j = 0; j < m_sizeJ; j++)
        {
            Vertex samplePosition = Vertex(i, 0.5f + j);
            float advectedVelocity = m_advectedVelocity.getU(i,j);
            float particleVelocity = m_fluidVelocityGrid.getU(i,j);
            float sdf = m_fluidSdf.lerpolateAt(samplePosition);
            m_fluidVelocityGrid.setU(i,j,nbcombine(advectedVelocity, particleVelocity, sdf, samplePosition));
        }
    }

    for (int i = 0; i < m_sizeI; i++)
    {
        for (int j = 0; j < m_sizeJ + 1; j++)
        {
            Vertex samplePosition = Vertex(0.5f + i, j);
            float advectedVelocity = m_advectedVelocity.getV(i,j);
            float particleVelocity = m_fluidVelocityGrid.getV(i,j);
            float sdf = m_fluidSdf.lerpolateAt(samplePosition);
            m_fluidVelocityGrid.setV(i,j,nbcombine(advectedVelocity, particleVelocity, sdf, samplePosition));
        }
    }
}

const Grid2d<float> FlipSmokeSolver::smokeConcentration() const
{
    return m_smokeConcentration;
}

const Grid2d<float> FlipSmokeSolver::temperature() const
{
    return m_temperature;
}

void FlipSmokeSolver::initAdditionalParameters()
{
    m_testValuePropertyIndex = m_markerParticles.addParticleProperty<float>();
    m_concentrationIndex = m_markerParticles.addParticleProperty<float>();
    m_temperatureIndex = m_markerParticles.addParticleProperty<float>();
}
