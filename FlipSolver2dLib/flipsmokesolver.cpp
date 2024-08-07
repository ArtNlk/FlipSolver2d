#include "flipsmokesolver.h"

#include "flipsolver2d.h"
#include "grid2d.h"
#include "mathfuncs.h"

#include <cmath>

FlipSmokeSolver::FlipSmokeSolver(const SmokeSolverParameters* p):
    FlipSolver(p),
    m_temperature(p->gridSizeI, p->gridSizeJ, p->ambientTemperature, OOBStrategy::OOB_CONST, p->ambientTemperature),
    m_smokeConcentration(p->gridSizeI, p->gridSizeJ, 0.f, OOBStrategy::OOB_CONST, 0.f),
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
    for (ssize_t i = 0; i < m_sizeI + 1; i++)
    {
        for (ssize_t j = 0; j < m_sizeJ + 1; j++)
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

    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_viscosityGrid.linearSize(),1,8);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSmokeSolver::centeredParamsToGridThread,this,range,
                                 std::ref(centeredWeights));
    }
    ThreadPool::i()->wait();
}

void FlipSmokeSolver::calcPressureRhs(std::vector<double> &rhs)
{
    const double scale = 1.f/m_dx;

    for (ssize_t i = 0; i < m_sizeI; i++)
    {
        for (ssize_t j = 0; j < m_sizeJ; j++)
        {
            if (!m_materialGrid.isSolid(i,j))
            {
                double val = -scale * divergenceAt(i,j);

                val -= m_materialGrid.isSolid(i-1,j) *
                       scale * static_cast<double>(m_fluidVelocityGrid.u(i,j) - 0);
                val += m_materialGrid.isSolid(i+1,j) *
                       scale * static_cast<double>(m_fluidVelocityGrid.u(i+1,j) - 0);
                val -= m_materialGrid.isSolid(i,j-1) *
                       scale * static_cast<double>(m_fluidVelocityGrid.v(i,j) - 0);
                val += m_materialGrid.isSolid(i,j+1) *
                       scale * static_cast<double>(m_fluidVelocityGrid.v(i,j+1) - 0);

                rhs.at(linearIndex(i,j)) = val;
            }
            else
            {
                rhs.at(linearIndex(i,j)) = 0.0;
            }
        }
    }
}

void FlipSmokeSolver::particleUpdate()
{
    FlipSolver::particleUpdate();


    for(size_t binIdx = 0; binIdx < m_markerParticles.bins().linearSize(); binIdx++)
    {
        ParticleBin& currentBin = m_markerParticles.bins().data().at(binIdx);
        std::vector<float>& particleTemperatures = currentBin.particleProperties<float>
                                                   (m_temperatureIndex);

        std::vector<float>& particleConcentrations = currentBin.particleProperties<float>
                                                     (m_concentrationIndex);

        for(size_t particleIdx = 0; particleIdx < currentBin.size(); particleIdx++)
        {
            //Vertex &position = currentBin.particlePosition(particleIdx);
            //particleTemperatures.at(i) = m_temperature.lerpolateAt(position);
            particleTemperatures.at(particleIdx) = m_ambientTemperature +
                    (particleTemperatures.at(particleIdx) - m_ambientTemperature) *
                    std::exp(-m_temperatureDecayRate * m_stepDt);
            //particleConcentrations.at(i) = m_smokeConcentration.lerpolateAt(position);
            particleConcentrations.at(particleIdx) *= std::exp(-m_concentrationDecayRate * m_stepDt);
            //p.testValue = p.smokeConcentrartion;
        }
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

    // std::vector<float>& particleViscosities = m_markerParticles.particleProperties<float>
    //                                           (m_viscosityPropertyIndex);

    for (ssize_t i = 0; i < m_sizeI; i++)
    {
        for (ssize_t j = 0; j < m_sizeJ; j++)
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

            ParticleBin& bin = m_markerParticles.binForGridIdx(i,j);
            auto& tempVec = bin.particleProperties<float>(m_temperatureIndex);
            auto& concVec = bin.particleProperties<float>(m_concentrationIndex);

            if(m_materialGrid.isSource(i,j))
            {
                int emitterId = m_emitterId.at(i,j);
                for(int p = 0; p < additionalParticles; p++)
                {
                    Vertex pos = jitteredPosInCell(i,j);
                    Vertex velocity = Vertex();
                    if(m_sources[emitterId].velocityTransfer())
                    {
                        velocity = m_fluidVelocityGrid.velocityAt(pos);
                    }
                    //                    int emitterId = m_emitterId.at(i,j);
                    //                    //float viscosity = m_sources[emitterId].viscosity();
                    //                    float conc = m_sources[emitterId].concentrartion();
                    //                    float temp = m_sources[emitterId].temperature();
                    float conc = m_smokeConcentration.interpolateAt(pos);
                    float temp = m_temperature.interpolateAt(pos);
                    size_t pIdx = bin.addMarkerParticle(pos,velocity);
                    tempVec[pIdx] = temp;
                    concVec[pIdx] = conc;
                    //testValues[pIdx] = temp / 1000.f;
                }
            }
        }
    }
}

void FlipSmokeSolver::gridUpdate()
{
    particleToGrid();
    m_stats.endStage(PARTICLE_TO_GRID);
    updateSdf();
    updateMaterials();
}

void FlipSmokeSolver::eulerAdvectParameters()
{
    Vertex offsetCentered(0.5f,0.5f);
    Grid2d<float> advectedConcentration(m_sizeI, m_sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetCentered);
    Grid2d<float> advectedTemperature(m_sizeI, m_sizeJ, 0.f, OOBStrategy::OOB_EXTEND, 0.f, offsetCentered);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(advectedConcentration.data().size());
    for(Range r : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSolver::eulerAdvectionThread,this,r,
                                 std::ref(m_smokeConcentration),std::ref(advectedConcentration));
        ThreadPool::i()->enqueue(&FlipSolver::eulerAdvectionThread,this,r,
                                 std::ref(m_temperature),std::ref(advectedTemperature));
    }
    ThreadPool::i()->wait();

    m_smokeConcentration = advectedConcentration;
    m_temperature = advectedTemperature;
    // for(int i = 0; i < m_smokeConcentration.linearSize(); i++)
    // {
    //     m_testGrid.data().at(i) = m_smokeConcentration.data().at(i);
    // }
}

void FlipSmokeSolver::applyPressuresToVelocityField(const std::vector<double> &pressures)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(pressures.size());

    for(Range range : ranges)
    {
        ThreadPool::i()->enqueue(&FlipSmokeSolver::applyPressureThreadU,this,range,std::ref(pressures));
        ThreadPool::i()->enqueue(&FlipSmokeSolver::applyPressureThreadV,this,range,std::ref(pressures));
    }
    ThreadPool::i()->wait();

    for(size_t i = 0; i < pressures.size(); i++)
    {
        m_testGrid.data().at(i) = pressures.at(i) / 100.0;
    }

//     if(anyNanInf(m_fluidVelocityGrid.velocityGridU().data()))
//     {
//         std::cout << "NaN or inf in U vector!\n" << std::flush;
//     }

//     if(anyNanInf(m_fluidVelocityGrid.velocityGridV().data()))
//     {
//         std::cout << "NaN or inf in V vector!\n" << std::flush;
//     }
}

void FlipSmokeSolver::applyPressureThreadU(Range range, const std::vector<double> &pressures)
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx);

    for (size_t pressureIdx = range.start; pressureIdx < range.end; pressureIdx++)
    {
        Index2d i2d = index2d(pressureIdx);
        Index2d i2dim1 = Index2d(i2d.i - 1,   i2d.j);
        //Index2d i2djm1 = Index2d(i2d.m_i,       i2d.m_j - 1);
        ssize_t pressureIndexIm1 = linearIndex(i2dim1);
        double pCurrent = pressures.at(pressureIdx);
        double pIm1 = pressureIndexIm1 == -1 ? 0.0 : pressures.at(pressureIndexIm1);
        //U part
        if(!m_materialGrid.isSolid(i2dim1) || !m_materialGrid.isSolid(i2d))
        {
            if(m_materialGrid.isSolid(i2dim1) || m_materialGrid.isSolid(i2d))
            {
                m_fluidVelocityGrid.setU(i2d,0.f);//Solids are stationary
            }
            else
            {
                m_fluidVelocityGrid.u(i2d) -= scale * (pCurrent - pIm1);
            }
        }
        else
        {
            m_fluidVelocityGrid.setUValidity(i2d,false);
        }
    }
}

void FlipSmokeSolver::applyPressureThreadV(Range range, const std::vector<double> &pressures)
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx);

    for (size_t pressureIdx = range.start; pressureIdx < range.end; pressureIdx++)
    {
        Index2d i2d = index2d(pressureIdx);
        Index2d i2djm1 = Index2d(i2d.i,   i2d.j - 1);
        //Index2d i2djm1 = Index2d(i2d.m_i,       i2d.m_j - 1);
        ssize_t pressureIndexJm1 = linearIndex(i2djm1);
        double pCurrent = pressures.at(pressureIdx);
        double pJm1 = pressureIndexJm1 == -1 ? 0.0 : pressures.at(pressureIndexJm1);
        //U part
        if(!m_materialGrid.isSolid(i2djm1) || !m_materialGrid.isSolid(i2d))
        {
            if(m_materialGrid.isSolid(i2djm1) || m_materialGrid.isSolid(i2d))
            {
                m_fluidVelocityGrid.setV(i2d,0.f);//Solids are stationary
            }
            else
            {
                m_fluidVelocityGrid.v(i2d) -= scale * (pCurrent - pJm1);
            }
        }
        else
        {
            m_fluidVelocityGrid.setVValidity(i2d,false);
        }
    }
}

void FlipSmokeSolver::seedInitialFluid()
{
    // std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float> (m_concentrationIndex);
    // std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>(m_temperatureIndex);

    for (ssize_t i = 0; i < m_sizeI; i++)
    {
        for (ssize_t j = 0; j < m_sizeJ; j++)
        {
            if (m_materialGrid.isStrictFluid(i, j))
            {
                ParticleBin& bin = m_markerParticles.binForGridIdx(i,j);
                auto& tempVec = bin.particleProperties<float>(m_temperatureIndex);
                auto& concVec = bin.particleProperties<float>(m_concentrationIndex);

                for (int p = 0; p < m_particlesPerCell; p++)
                {
                    Vertex pos = jitteredPosInCell(i, j);
                    Vertex velocity = m_fluidVelocityGrid.velocityAt(pos);
                    float conc = m_smokeConcentration.interpolateAt(pos);
                    float temp = m_temperature.interpolateAt(pos);
                    size_t pIdx = bin.addMarkerParticle(pos,velocity);
                    tempVec[pIdx] = temp;
                    concVec[pIdx] = conc;
                }
            }
        }
    }
}

IndexedPressureParameters FlipSmokeSolver::getPressureProjectionMatrix()
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);

    IndexedPressureParameters output(linearSize()*0.33, *this, scale);

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    std::vector<Range> threadRanges = ThreadPool::i()->splitRange(linearSize());
    size_t currRangeIdx = 0;

    for(ssize_t i = 0; i < m_sizeI; i++)
    {
        for(ssize_t j = 0; j < m_sizeJ; j++)
        {
            const ssize_t linIdx = indexer.linearIndex(i,j);

            if(m_materialGrid.isSolid(i,j))
            {
                if(linIdx >= threadRanges.at(currRangeIdx).end)
                {
                    output.endThreadDataRange();
                    currRangeIdx++;
                }
                continue;
            }

            IndexedPressureParameterUnit unit;
            unit.unitIndex = linIdx;

            const ssize_t linIdxPosAx = indexer.linearIdxOfOffset(linIdx,1,0);
            const ssize_t linIdxPosAy = indexer.linearIdxOfOffset(linIdx,0,1);
            const ssize_t linIdxNegAx = indexer.linearIdxOfOffset(linIdx,-1,0);
            const ssize_t linIdxNegAy = indexer.linearIdxOfOffset(linIdx,0,-1);

            double diag = 0.0;
            //X Neighbors
            if(m_materialGrid.isFluid(i-1,j))
            {
                unit.nonsolidNeighborCount++;
                unit.setNeighbor(I_NEG_NEIGHBOR_BIT, inBounds(linIdxNegAx));
            }else if(m_materialGrid.isEmpty(i-1,j))
            {
                unit.nonsolidNeighborCount++;
            }

            if(m_materialGrid.isFluid(i+1,j))
            {
                unit.nonsolidNeighborCount++;
                unit.setNeighbor(I_POS_NEIGHBOR_BIT, inBounds(linIdxPosAx));
            } else if(m_materialGrid.isEmpty(i+1,j))
            {
                unit.nonsolidNeighborCount++;
            }

            //Y Neighbors
            if(m_materialGrid.isFluid(i,j-1))
            {
                unit.nonsolidNeighborCount++;
                unit.setNeighbor(J_NEG_NEIGHBOR_BIT, inBounds(linIdxNegAy));
            }else if(m_materialGrid.isEmpty(i,j-1))
            {
                unit.nonsolidNeighborCount++;
            }

            if(m_materialGrid.isFluid(i,j+1))
            {
                unit.nonsolidNeighborCount++;
                unit.setNeighbor(J_POS_NEIGHBOR_BIT, inBounds(linIdxPosAy));
            } else if(m_materialGrid.isEmpty(i,j+1))
            {
                unit.nonsolidNeighborCount++;
            }

            if(linIdx >= threadRanges.at(currRangeIdx).end)
            {
                output.endThreadDataRange();
                currRangeIdx++;
            }

            output.add(unit);
        }
    }

    if(currRangeIdx != threadRanges.size())
    {
        output.endThreadDataRange();
    }

    return output;
}

IndexedIPPCoefficients FlipSmokeSolver::getIPPCoefficients(const IndexedPressureParameters &mat)
{
    const double scale = m_stepDt / (m_fluidDensity * m_dx * m_dx);

    IndexedIPPCoefficients output(linearSize()*0.33, *this);

    LinearIndexable2d& indexer = *dynamic_cast<LinearIndexable2d*>(this);

    std::vector<Range> threadRanges = ThreadPool::i()->splitRange(linearSize());
    size_t currRangeIdx = 0;

    size_t matEntryIdx = 0;

    for(ssize_t i = 0; i < m_sizeI; i++)
    {
        for(ssize_t j = 0; j < m_sizeJ; j++)
        {
            const ssize_t linIdx = indexer.linearIndex(i,j);

            if(m_materialGrid.isSolid(i,j))
            {
                if(linIdx >= threadRanges.at(currRangeIdx).end)
                {
                    output.endThreadDataRange();
                    currRangeIdx++;
                }
                continue;
            }

            IndexedIPPCoefficientUnit unit;
            unit.unitIndex = linIdx;

            const double invscale = scale/(mat.data().at(matEntryIdx).nonsolidNeighborCount*scale);

            const ssize_t jNegLinIdx = indexer.linearIdxOfOffset(linIdx,0,-1);
            const ssize_t jNegP1LinIdx = indexer.linearIdxOfOffset(linIdx,0,-1)+1;
            const ssize_t iNegLinIdx = indexer.linearIdxOfOffset(linIdx,-1,0);
            const ssize_t iPosLinIdx = indexer.linearIdxOfOffset(linIdx,1,0);
            const ssize_t jPosM1LinIdx = indexer.linearIdxOfOffset(linIdx,0,1)-1;
            const ssize_t jPosLinIdx = indexer.linearIdxOfOffset(linIdx,0,1);

            if(jNegLinIdx >= 0 && m_materialGrid.isFluid(jNegLinIdx))
            {
                unit.jNeg = invscale;
            }

            if(jNegP1LinIdx >= 0 && m_materialGrid.isFluid(jNegP1LinIdx))
            {
                unit.jNegP1 = invscale*invscale;
            }

            if(iNegLinIdx >= 0 && m_materialGrid.isFluid(iNegLinIdx))
            {
                unit.iNeg = invscale;
            }

            unit.diag = (invscale*invscale)*2.0 + 1.0;

            if(iPosLinIdx >= 0 && m_materialGrid.isFluid(iPosLinIdx))
            {
                unit.iPos = invscale;
            }

            if(jPosM1LinIdx >= 0 && m_materialGrid.isFluid(jPosM1LinIdx))
            {
                unit.jPosM1 = invscale*invscale;
            }

            if(jPosLinIdx >= 0 && m_materialGrid.isFluid(jPosLinIdx))
            {
                unit.jPos = invscale;
            }

            if(linIdx >= threadRanges.at(currRangeIdx).end)
            {
                output.endThreadDataRange();
                currRangeIdx++;
            }

            output.add(unit);
            matEntryIdx++;
        }
    }

    if(currRangeIdx != threadRanges.size())
    {
        output.endThreadDataRange();
    }

    return output;
}

void FlipSmokeSolver::centeredParamsToGridThread(Range r, Grid2d<float> &cWeights)
{
    // std::vector<float>& particleTemperatures = m_markerParticles.particleProperties<float>
    //                                            (m_temperatureIndex);

    // std::vector<float>& particleConcentrations = m_markerParticles.particleProperties<float>
    //                                              (m_concentrationIndex);

    for(size_t idx = r.start; idx < r.end; idx++)
    {
        Index2d i2d = m_fluidVelocityGrid.index2d(idx);
        std::array<ssize_t,9> affectingBins = m_markerParticles.binsForGridCell(i2d);
        for(ssize_t binIdx : affectingBins)
        {
            if(binIdx >= 0)
            {
                ParticleBin& bin = m_markerParticles.binForBinIdx(binIdx);
                auto& tempVec = bin.particleProperties<float>(m_temperatureIndex);
                auto& concVec = bin.particleProperties<float>(m_concentrationIndex);

                for(size_t particleIdx = 0; particleIdx < bin.size(); particleIdx++)
                {
                    Vertex position = bin.particlePosition(particleIdx);
                    float weightCentered = simmath::quadraticBSpline(position.x() - (i2d.i),
                                                                     position.y() - (i2d.j));
                    if(std::abs(weightCentered) > 1e-6f)
                    {
                        cWeights.at(i2d.i,i2d.j) += weightCentered;
                        m_temperature.at(i2d.i,i2d.j) +=
                            weightCentered * tempVec.at(particleIdx);

                        m_smokeConcentration.at(i2d.i,i2d.j) +=
                            weightCentered * concVec.at(particleIdx);
                        m_knownCenteredParams.at(i2d.i,i2d.j) = true;
                    }
                }
            }
        }

        if(m_knownCenteredParams.inBounds(i2d))
        {
            if(m_knownCenteredParams.at(i2d))
            {
                m_temperature.at(i2d) /= cWeights.at(i2d);
                m_smokeConcentration.at(i2d) /= cWeights.at(i2d);
                //m_testGrid.at(i,j) = m_viscosityGrid.at(i,j);
            }
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
    m_concentrationIndex = m_markerParticles.addParticleProperty<float>();
    m_temperatureIndex = m_markerParticles.addParticleProperty<float>();
}
