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
    m_particlePositions(0),
    m_velocities(0),
    m_markedForDeath(0),
    m_markedForRebin(0),
    m_rebinningSet(0),
    m_properties(0)
{

}

MarkerParticleSystem::ParticleBin &MarkerParticleSystem::binForGridIdx(int linIdx)
{
    return m_particleBins.data()[gridToBinIdx(linIdx)];
}

MarkerParticleSystem::ParticleBin &MarkerParticleSystem::binForGridIdx(Index2d idx)
{
    return binForGridIdx(idx.i, idx.j);
}

MarkerParticleSystem::ParticleBin &MarkerParticleSystem::binForGridIdx(int i, int j)
{
    return m_particleBins.data()[gridToBinIdx(i,j)];
}

MarkerParticleSystem::ParticleBin &MarkerParticleSystem::binForGridPosition(Vertex pos)
{
    return binForGridIdx(pos.x(), pos.y());
}

MarkerParticleSystem::ParticleBin &MarkerParticleSystem::binForBinIdx(int linIdx)
{
    return m_particleBins.data()[linIdx];
}

void MarkerParticleSystem::rebinParticles()
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(m_particleBins.linearSize());
    std::latch sync(ThreadPool::i()->threadCount());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&MarkerParticleSystem::rebinParticlesThread,this,range,std::ref(sync));
    }
    ThreadPool::i()->wait();
    m_rebinningSet.clear();
}

void MarkerParticleSystem::pruneParticles()
{
    for(int pIndex = 0; pIndex < m_particlePositions.size(); pIndex++)
    {
        if(m_markedForDeath[pIndex])
        {
            eraseMarkerParticle(pIndex);
            pIndex--;
        }
    }
}

size_t MarkerParticleSystem::addMarkerParticle(Vertex position, Vertex velocity)
{
    m_particlePositions.push_back(position);
    binForGridPosition(m_particlePositions.back()).push_back(m_particlePositions.size() - 1);
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

void MarkerParticleSystem::eraseMarkerParticle(size_t index)
{
    for(size_t rebIdx : m_rebinningSet)
    {
        if(rebIdx == index)
        {
            std::cout << "Attempt to delete rebinned particle!" << std::endl;
        }
    }
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

void MarkerParticleSystem::markForDeath(size_t particleIndex)
{
    m_markedForDeath.at(particleIndex) = true;
}

bool MarkerParticleSystem::markedForDeath(size_t particleIdx)
{
    return m_markedForDeath.at(particleIdx);
}

void MarkerParticleSystem::scheduleRebin(size_t particleIdx)
{
    m_rebinningSet.push_back(particleIdx);
    m_markedForRebin.at(particleIdx) = true;
}

bool MarkerParticleSystem::scheduledForRebin(size_t particleIdx)
{
    return m_markedForRebin[particleIdx];
}

Grid2d<MarkerParticleSystem::ParticleBin>& MarkerParticleSystem::bins()
{
    return m_particleBins;
}

Vertex &MarkerParticleSystem::particlePosition(size_t index)
{
    return m_particlePositions.at(index);
}

Vertex &MarkerParticleSystem::particleVelocity(size_t index)
{
    return m_velocities.at(index);
}

std::vector<Vertex> &MarkerParticleSystem::positions()
{
    return m_particlePositions;
}

std::vector<Vertex> &MarkerParticleSystem::velocities()
{
    return m_velocities;
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

std::array<int, 9> MarkerParticleSystem::rebinSetForBinIdx(Index2d idx)
{
    std::array<int,9> output({-1,-1,-1,-1,-1,-1,-1,-1,-1});
    int outputIdx = 0;
    for(int iOffset = -1; iOffset <= 1; iOffset++)
    {
        for(int jOffset = -1; jOffset <= 1; jOffset++)
        {
            if(m_particleBins.inBounds(idx.i + iOffset, idx.j + jOffset))
            {
                output[outputIdx] = m_particleBins.linearIndex(idx.i + iOffset,
                                                               idx.j + jOffset);
            }
            outputIdx++;
        }
    }

    return output;
}

size_t MarkerParticleSystem::particleCount() const
{
    return m_particlePositions.size();
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

void MarkerParticleSystem::rebinParticlesThread(Range r, std::latch &sync)
{
    std::vector<float>& testValues = particleProperties<float>(0);

    for(int binIdx = r.start; binIdx < r.end; binIdx++)
    {
        ParticleBin& currentBin = m_particleBins.data()[binIdx];
        for(size_t pIdx = 0; pIdx < currentBin.size(); pIdx++)
        {
            if(m_markedForRebin[currentBin[pIdx]])
            {
                m_markedForRebin[currentBin[pIdx]] = false;
                swapErase(currentBin,pIdx);
                pIdx--;
            }
        }
        //currentBin.clear();
        for(const size_t particleIndex : m_rebinningSet)
        {
            int idx = gridToBinIdx(m_particlePositions[particleIndex]);

            if(idx == binIdx)
            {
                currentBin.push_back(particleIndex);
            }
        }
    }
}


