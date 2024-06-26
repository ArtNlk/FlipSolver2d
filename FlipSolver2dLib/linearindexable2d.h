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
    LinearIndexable2d(int sizeI, int sizeJ) :
        m_sizeI(sizeI),
        m_sizeJ(sizeJ)
    {
    }

    inline int sizeI() const
    {
        return m_sizeI;
    }

    inline int sizeJ() const
    {
        return m_sizeJ;
    }

    inline int linearIndex(int i, int j) const
    {
        if(i < 0 || i >= m_sizeI || j < 0 || j >= m_sizeJ)
        {
            return -1;
        }
        return i * m_sizeJ + j;
    }

    inline int linearIndex(Index2d index) const
    {
        if(index.i < 0 || index.i >= m_sizeI || index.j < 0 || index.j >= m_sizeJ) return -1;
        return index.i * m_sizeJ + index.j;
    }

    inline Index2d index2d(int linearIndex) const
    {
        ASSERT_BETWEEN(linearIndex,-1,m_sizeJ*m_sizeI);
        Index2d index;
        index.i = linearIndex / m_sizeJ;
        index.j = linearIndex - index.i*m_sizeJ;

        return index;
    }

    inline std::array<int,8> getNeighborhood(int index)
    {
        return getNeighborhood(index2d(index));
    }

    inline std::array<int, 8> getNeighborhood(Index2d index)
    {
        return getNeighborhood(index.i, index.j);
    }

    inline std::array<int, 8> getNeighborhood(int i, int j)
    {
        std::array<int, 8> output;
        int outputIdx = 0;
        for(int iOffset = -1; iOffset <= 1; iOffset++)
        {
            for(int jOffset = -1; jOffset <= 1; jOffset++)
            {
                if(iOffset == 0 && jOffset == 0) continue;
                int index = linearIndex(i + iOffset,j + jOffset);
                output[outputIdx] = index;
                outputIdx++;
            }
        }

        return output;
    }

    inline std::array<int,4> immidiateNeighbors(Index2d idx) const
    {
        return immidiateNeighbors(linearIndex(idx));
    }

    inline std::array<int,4> immidiateNeighbors(int i, int j) const
    {
        return immidiateNeighbors(linearIndex(i,j));
    }

    inline std::array<int,4> immidiateNeighbors(int linearIdx) const
    {
        std::array<int,4> output;
        output[0] = std::clamp(linearIdxOfOffset(linearIdx,-1,0),0,m_sizeI*m_sizeJ - 1);
        output[1] = std::clamp(linearIdxOfOffset(linearIdx,1,0),0,m_sizeI*m_sizeJ - 1);
        output[2] = std::clamp(linearIdxOfOffset(linearIdx,0,-1),0,m_sizeI*m_sizeJ - 1);
        output[3] = std::clamp(linearIdxOfOffset(linearIdx,0,1),0,m_sizeI*m_sizeJ - 1);
        return output;
    }

    inline int linearIdxOfOffset(int linearIdx, int iOffset, int jOffset) const
    {
        return linearIdx + iOffset * m_sizeJ + jOffset;
    }

    inline size_t linearIdxOfOffset(size_t linearIdx, int iOffset, int jOffset) const
    {
        return linearIdx + iOffset * m_sizeJ + jOffset;
    }

    inline size_t iLinearOffset()
    {
        return m_sizeJ;
    }

    inline size_t jLinearOffset()
    {
        return 1;
    }

    inline bool inBounds(Index2d index) const
    {
        return (index.i >= 0 && index.i < m_sizeI && index.j >= 0 && index.j < m_sizeJ);
    }

    inline bool inBounds(int i, int j) const
    {
        return (i >= 0 && i < m_sizeI && j >= 0 && j < m_sizeJ);
    }

    inline bool inBounds(int linearIndex) const
    {
        return linearIndex >= 0 && linearIndex < m_sizeI * m_sizeJ;
    }

    size_t linearSize() const
    {
        return m_sizeI * m_sizeJ;
    }

protected:
    int m_sizeI;
    int m_sizeJ;
};

#endif // LINEARINDEXABLE2D_H
