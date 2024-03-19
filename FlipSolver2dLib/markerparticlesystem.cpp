#include "markerparticlesystem.h"
#include "threadpool.h"
#include <algorithm>
#include <cstddef>
#include <latch>
#include <variant>
#include <vector>

MarkerParticleSystem::MarkerParticleSystem(int gridSizeI, int gridSizeJ, size_t binSize):
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

ParticleBin &MarkerParticleSystem::binForGridIdx(int linIdx)
{
    return m_particleBins.data()[gridToBinIdx(linIdx)];
}

ParticleBin &MarkerParticleSystem::binForGridIdx(Index2d idx)
{
    return binForGridIdx(idx.i, idx.j);
}

ParticleBin &MarkerParticleSystem::binForGridIdx(int i, int j)
{
    return m_particleBins.data()[gridToBinIdx(i,j)];
}

ParticleBin &MarkerParticleSystem::binForGridPosition(Vertex pos)
{
    return binForGridIdx(pos.x(), pos.y());
}

ParticleBin &MarkerParticleSystem::binForBinIdx(int linIdx)
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

void MarkerParticleSystem::scheduleRebin(size_t binIdx, size_t particleIdx)
{
    m_particleBins.data().at(binIdx).scheduleRebin(m_rebinningSet,particleIdx);
}

void MarkerParticleSystem::pruneParticles()
{
    for(ParticleBin& bin : m_particleBins.data())
    {
        bin.pruneParticles();
    }
}

Grid2d<ParticleBin>& MarkerParticleSystem::bins()
{
    return m_particleBins;
}

void MarkerParticleSystem::addMarkerParticle(size_t binIdx, Vertex position, Vertex velocity)
{
    ASSERT(binIdx < m_bins.linearSize());

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

Index2d MarkerParticleSystem::binIdxForIdx(int i, int j)
{
    int binI = i/m_binSize;
    int binJ = j/m_binSize;
    return Index2d(binI, binJ);
}

std::array<int, 9> MarkerParticleSystem::binsForGridCell(Index2d idx)
{
    return binsForGridCell(idx.i, idx.j);
}

std::array<int, 9> MarkerParticleSystem::binsForGridCell(int i, int j)
{
    std::array<int,9> output({-1,-1,-1,-1,-1,-1,-1,-1,-1});
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

int MarkerParticleSystem::gridToBinIdx(Index2d idx)
{
    return gridToBinIdx(idx.i,idx.j);
}

int MarkerParticleSystem::gridToBinIdx(Vertex pos)
{
    return gridToBinIdx(pos.x(), pos.y());
}

int MarkerParticleSystem::gridToBinIdx(int linIdx)
{
    return gridToBinIdx(m_gridIndexer.index2d(linIdx));
}

int MarkerParticleSystem::gridToBinIdx(int i, int j)
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
    return m_particlePositions.at(index);
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

void ParticleBin::rebinParticles(ParticleBin &rebinningSet)
{
    for(size_t i = 0; i < rebinningSet.size(); i++)
    {
        if(rebinningSet.binIdx() == m_binIdx)
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

void ParticleBin::scheduleRebin(ParticleBin &rebinningSet, size_t idx)
{
    moveParticleToBin(rebinningSet, idx);
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
