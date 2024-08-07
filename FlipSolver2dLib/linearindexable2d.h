#ifndef LINEARINDEXABLE2D_H
#define LINEARINDEXABLE2D_H

#include <algorithm>
#include <array>
#include <vector>

#include "index2d.h"
#include "customassert.h"

class LinearIndexable2d
{
public:
    LinearIndexable2d(size_t sizeI, size_t sizeJ) :
        m_sizeI(sizeI),
        m_sizeJ(sizeJ)
    {
    }

    ssize_t sizeI() const
    {
        return m_sizeI;
    }

    ssize_t sizeJ() const
    {
        return m_sizeJ;
    }

    ssize_t linearIndex(ssize_t i, ssize_t j) const
    {
        if(i < 0 || i >= m_sizeI || j < 0 || j >= m_sizeJ)
        {
            return -1;
        }
        return i * m_sizeJ + j;
    }

    ssize_t linearIndex(Index2d index) const
    {
        if(index.i < 0 || index.i >= m_sizeI || index.j < 0 || index.j >= m_sizeJ) return -1;
        return index.i * m_sizeJ + index.j;
    }

    Index2d index2d(ssize_t linearIndex) const
    {
        ASSERT_BETWEEN(linearIndex,-1,m_sizeJ*m_sizeI);
        Index2d index;
        index.i = linearIndex / m_sizeJ;
        index.j = linearIndex - index.i*m_sizeJ;

        return index;
    }

    std::array<ssize_t,8> getNeighborhood(size_t index)
    {
        return getNeighborhood(index2d(index));
    }

    std::array<ssize_t, 8> getNeighborhood(Index2d index)
    {
        return getNeighborhood(index.i, index.j);
    }

    std::array<ssize_t, 8> getNeighborhood(size_t i, size_t j)
    {
        std::array<ssize_t, 8> output;
        size_t outputIdx = 0;
        for(ssize_t iOffset = -1; iOffset <= 1; iOffset++)
        {
            for(ssize_t jOffset = -1; jOffset <= 1; jOffset++)
            {
                if(iOffset == 0 && jOffset == 0) continue;
                ssize_t index = linearIndex(i + iOffset,j + jOffset);
                output[outputIdx] = index;
                outputIdx++;
            }
        }

        return output;
    }

    std::array<ssize_t,4> immidiateNeighbors(Index2d idx) const
    {
        return immidiateNeighbors(linearIndex(idx));
    }

    std::array<ssize_t,4> immidiateNeighbors(ssize_t i, ssize_t j) const
    {
        return immidiateNeighbors(linearIndex(i,j));
    }

    std::array<ssize_t,4> immidiateNeighbors(ssize_t linearIdx) const
    {
        std::array<ssize_t,4> output;
        output[0] = std::clamp<ssize_t>(linearIdxOfOffset(linearIdx,-1,0),0,m_sizeI*m_sizeJ - 1);
        output[1] = std::clamp<ssize_t>(linearIdxOfOffset(linearIdx,1,0),0,m_sizeI*m_sizeJ - 1);
        output[2] = std::clamp<ssize_t>(linearIdxOfOffset(linearIdx,0,-1),0,m_sizeI*m_sizeJ - 1);
        output[3] = std::clamp<ssize_t>(linearIdxOfOffset(linearIdx,0,1),0,m_sizeI*m_sizeJ - 1);
        return output;
    }

    ssize_t linearIdxOfOffset(ssize_t linearIdx, ssize_t iOffset, ssize_t jOffset) const
    {
        return linearIdx + iOffset * m_sizeJ + jOffset;
    }

    size_t iLinearOffset()
    {
        return m_sizeJ;
    }

    size_t jLinearOffset()
    {
        return 1;
    }

    bool inBounds(Index2d index) const
    {
        return (index.i >= 0 && index.i < m_sizeI && index.j >= 0 && index.j < m_sizeJ);
    }

    bool inBounds(ssize_t i, ssize_t j) const
    {
        return (i >= 0 && i < m_sizeI && j >= 0 && j < m_sizeJ);
    }

    bool inBounds(ssize_t linearIndex) const
    {
        return linearIndex >= 0 && linearIndex < m_sizeI * m_sizeJ;
    }

    size_t linearSize() const
    {
        return m_sizeI * m_sizeJ;
    }

protected:
    ssize_t m_sizeI;
    ssize_t m_sizeJ;
};

#endif // LINEARINDEXABLE2D_H
