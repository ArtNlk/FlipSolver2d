#include "markerparticlesystem.h"
#include "threadpool.h"
#include <algorithm>
#include <cstddef>
#include <latch>
#include <variant>
#include <vector>

MarkerParticleSystem::MarkerParticleSystem(size_t gridSizeI, size_t gridSizeJ, size_t binSize):
    m_binSize(binSize),
    m_gridIndexer(gridSizeI, gridSizeJ),
    m_particleBins((gridSizeI / binSize) + (gridSizeI % binSize != 0),
                     (gridSizeJ / binSize) + (gridSizeJ % binSize != 0)),
    m_rebinningSet()
{
    for(size_t i = 0; i < m_particleBins.linearSize(); i++)
    {
        m_particleBins.data()[i].setBinIndex(i);
    }
}

ParticleBin &MarkerParticleSystem::binForGridIdx(size_t linIdx)
{
    return m_particleBins.data()[gridToBinIdx(linIdx)];
}

ParticleBin &MarkerParticleSystem::binForGridIdx(Index2d idx)
{
    return binForGridIdx(idx.i, idx.j);
}

ParticleBin &MarkerParticleSystem::binForGridIdx(size_t i, size_t j)
{
    return m_particleBins.data()[gridToBinIdx(i,j)];
}

ParticleBin &MarkerParticleSystem::binForGridPosition(Vertex pos)
{
    return binForGridIdx(pos.x(), pos.y());
}

ParticleBin &MarkerParticleSystem::binForBinIdx(size_t linIdx)
{
    return m_particleBins.data()[linIdx];
}

void MarkerParticleSystem::rebinParticles()
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_particleBins.linearSize());
    std::latch sync(ThreadPool::i()->threadCount());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&MarkerParticleSystem::rebinParticlesThread,this,range);
    }
    ThreadPool::i()->wait();
    m_rebinningSet.clear();
}

void MarkerParticleSystem::scheduleRebin(size_t binIdx, RebinRecord r)
{
    m_particleBins.data().at(binIdx).scheduleRebin(m_rebinningSet,r.particleIdx, r.newBinIdx);
}

Grid2d<ParticleBin>& MarkerParticleSystem::bins()
{
    return m_particleBins;
}

void MarkerParticleSystem::addMarkerParticle(size_t binIdx, Vertex position, Vertex velocity)
{
    ASSERT(binIdx < m_particleBins.linearSize());

    m_particleBins.data().at(binIdx).addMarkerParticle(position, velocity);
}

void MarkerParticleSystem::markForDeath(size_t binIdx, size_t particleIdx)
{
    m_particleBins.data().at(binIdx).markForDeath(particleIdx);
}

std::vector<Vertex> &MarkerParticleSystem::positions(size_t binIdx)
{
    return m_particleBins.data().at(binIdx).positions();
}

std::vector<Vertex> &MarkerParticleSystem::velocities(size_t binIdx)
{
    return m_particleBins.data().at(binIdx).velocities();
}

Index2d MarkerParticleSystem::binIdxForIdx(Index2d idx)
{
    return binIdxForIdx(idx.i, idx.j);
}

Index2d MarkerParticleSystem::binIdxForIdx(ssize_t i, ssize_t j)
{
    int binI = i/m_binSize;
    int binJ = j/m_binSize;
    return Index2d(binI, binJ);
}

std::array<ssize_t, 9> MarkerParticleSystem::binsForGridCell(Index2d idx)
{
    return binsForGridCell(idx.i, idx.j);
}

std::array<ssize_t, 9> MarkerParticleSystem::binsForGridCell(ssize_t i, ssize_t j)
{
    std::array<ssize_t,9> output({-1,-1,-1,-1,-1,-1,-1,-1,-1});
    int outputIdx = 0;
    Index2d centerBinIdx = binIdxForIdx(i,j);
    for(int iOffset = -1; iOffset <= 1; iOffset++)
    {
        for(int jOffset = -1; jOffset <= 1; jOffset++)
        {
            if(m_particleBins.inBounds(centerBinIdx.i + iOffset, centerBinIdx.j + jOffset))
            {
                output[outputIdx] = m_particleBins.linearIndex(centerBinIdx.i + iOffset,
                                                               centerBinIdx.j + jOffset);
            }
            outputIdx++;
        }
    }

    return output;
}

size_t MarkerParticleSystem::particleCount() const
{
    size_t count = 0;
    for(const ParticleBin& bin : m_particleBins.data())
    {
        count += bin.size();
    }

    return count;
}

ssize_t MarkerParticleSystem::gridToBinIdx(Index2d idx)
{
    return gridToBinIdx(idx.i,idx.j);
}

ssize_t MarkerParticleSystem::gridToBinIdx(Vertex pos)
{
    return gridToBinIdx(pos.x(), pos.y());
}

ssize_t MarkerParticleSystem::gridToBinIdx(ssize_t linIdx)
{
    return gridToBinIdx(m_gridIndexer.index2d(linIdx));
}

ssize_t MarkerParticleSystem::gridToBinIdx(ssize_t i, ssize_t j)
{
    return m_particleBins.linearIndex(binIdxForIdx(i,j));
}

void MarkerParticleSystem::rebinParticlesThread(Range r)
{
    for(int i = r.start; i < r.end; i++)
    {
        m_particleBins.data().at(i).rebinParticles(m_rebinningSet);
    }
}

ParticleBin::ParticleBin(size_t binIdx) :
    m_binIdx(binIdx)
{

}

size_t ParticleBin::addMarkerParticle(Vertex position, Vertex velocity)
{
    m_particlePositions.push_back(position);
    m_velocities.push_back(velocity);
    m_markedForDeath.push_back(false);
    m_markedForRebin.push_back(false);

    for(VariantVector& v : m_properties)
    {
        switch(v.index())
        {
        case ParticleAttributeType::ATTR_FLOAT:
            std::get<std::vector<float>>(v).push_back(0.f);
            break;

        default:
            break;
        }
    }

    return m_particlePositions.size() - 1;
}

void ParticleBin::eraseMarkerParticle(size_t index)
{
    swapErase(m_particlePositions, index);
    swapErase(m_velocities, index);
    swapErase(m_markedForDeath, index);
    swapErase(m_markedForRebin, index);
    for(VariantVector& v : m_properties)
    {
        switch(v.index())
        {
        case ParticleAttributeType::ATTR_FLOAT:
            swapErase(std::get<std::vector<float>>(v),index);
            break;

        default:
            break;
        }
    }
}

void ParticleBin::markForDeath(size_t index)
{
    m_markedForDeath.at(index) = true;
}

bool ParticleBin::markedForDeath(size_t index)
{
    return m_markedForDeath.at(index);
}

Vertex &ParticleBin::particlePosition(size_t index)
{
    return m_particlePositions.at(index);
}

Vertex &ParticleBin::particleVelocity(size_t index)
{
    return m_velocities.at(index);
}

void ParticleBin::pruneParticles()
{
    for(size_t i = 0; i < m_particlePositions.size(); i++)
    {
        if(m_markedForDeath[i])
        {
            eraseMarkerParticle(i);
            i--;
        }
    }
}

void ParticleBin::rebinParticles(RebinSet &rebinningSet)
{
    for(size_t i = 0; i < rebinningSet.size(); i++)
    {
        if(rebinningSet.particleBinIdx(i) == m_binIdx)
        {
            rebinningSet.copyParticleToBin(*this, i);
        }
    }
}

void ParticleBin::moveParticleToBin(ParticleBin &other, size_t idx)
{
    ASSERT(other.m_properties.size() == m_properties.size());
    other.m_particlePositions.push_back(swapPop(m_particlePositions,idx));
    other.m_velocities.push_back(swapPop(m_velocities,idx));
    other.m_markedForDeath.push_back(swapPop(m_markedForDeath,idx));
    other.m_markedForRebin.push_back(swapPop(m_markedForRebin,idx));
    for(int i = 0; i < m_properties.size(); i++)
    {
        VariantVector& v = m_properties.at(i);
        VariantVector& otherV = other.m_properties.at(i);
        switch(v.index())
        {
        case ParticleAttributeType::ATTR_FLOAT:
            std::get<std::vector<float>>(otherV).push_back(swapPop(std::get<std::vector<float>>(v),idx));
            break;

        default:
            break;
        }
    }
}

void ParticleBin::copyParticleToBin(ParticleBin &other, size_t idx)
{
    ASSERT(other.m_properties.size() == m_properties.size());
    other.m_particlePositions.push_back(m_particlePositions.at(idx));
    other.m_velocities.push_back(m_velocities.at(idx));
    other.m_markedForDeath.push_back(m_markedForDeath.at(idx));
    other.m_markedForRebin.push_back(m_markedForRebin.at(idx));
    for(int i = 0; i < m_properties.size(); i++)
    {
        VariantVector& v = m_properties.at(i);
        VariantVector& otherV = other.m_properties.at(i);
        switch(v.index())
        {
        case ParticleAttributeType::ATTR_FLOAT:
            std::get<std::vector<float>>(otherV).push_back(std::get<std::vector<float>>(v).at(idx));
            break;

        default:
            break;
        }
    }
}

void ParticleBin::scheduleRebin(RebinSet &rebinningSet, size_t particleIndex, size_t newBinIndex)
{
    ASSERT(rebinningSet.m_properties.size() == m_properties.size());
    rebinningSet.m_particlePositions.push_back(m_particlePositions.at(particleIndex));
    rebinningSet.m_velocities.push_back(m_velocities.at(particleIndex));
    rebinningSet.m_markedForDeath.push_back(m_markedForDeath.at(particleIndex));
    rebinningSet.m_markedForRebin.push_back(m_markedForRebin.at(particleIndex));

    rebinningSet.binIndexes().push_back(newBinIndex);

    for(int i = 0; i < m_properties.size(); i++)
    {
        VariantVector& v = m_properties.at(i);
        VariantVector& otherV = rebinningSet.m_properties.at(i);
        switch(v.index())
        {
        case ParticleAttributeType::ATTR_FLOAT:
            std::get<std::vector<float>>(otherV).push_back(std::get<std::vector<float>>(v).at(particleIndex));
            break;

        default:
            break;
        }
    }
    markForDeath(particleIndex);
}

void ParticleBin::clear()
{
    m_particlePositions.clear();
    m_velocities.clear();
    m_markedForDeath.clear();
    m_markedForRebin.clear();

    for(int i = 0; i < m_properties.size(); i++)
    {
        VariantVector& v = m_properties.at(i);
        switch(v.index())
        {
        case ParticleAttributeType::ATTR_FLOAT:
            std::get<std::vector<float>>(v).clear();
            break;

        default:
            break;
        }
    }
}

void ParticleBin::setBinIndex(size_t index)
{
    m_binIdx = index;
}
