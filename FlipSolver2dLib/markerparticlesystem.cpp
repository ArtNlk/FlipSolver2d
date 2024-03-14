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
    m_rebinningSet(0)
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

}

void MarkerParticleSystem::addMarkerParticle(Vertex position, Vertex velocity)
{
    ParticleBin& bin = binForGridPosition(position);
    bin.emplace_back(position, velocity);
}

Grid2d<MarkerParticleSystem::ParticleBin>& MarkerParticleSystem::bins()
{
    return m_particleBins;
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

}


