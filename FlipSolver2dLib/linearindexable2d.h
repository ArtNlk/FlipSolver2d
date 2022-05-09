#ifndef LINEARINDEXABLE2D_H
#define LINEARINDEXABLE2D_H

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

public:
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
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        if(i < 0 || i >= m_sizeI || j < 0 || j >= m_sizeJ) return -1;
        return i * m_sizeJ + j;
    }

    inline int linearIndex(Index2d index) const
    {
        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
        if(index.m_i < 0 || index.m_i >= m_sizeI || index.m_j < 0 || index.m_j >= m_sizeJ) return -1;
        return index.m_i * m_sizeJ + index.m_j;
    }

    inline Index2d index2d(int linearIndex) const
    {
        ASSERT_BETWEEN(linearIndex,-1,m_sizeJ*m_sizeI);
        Index2d index;
        index.m_i = linearIndex / m_sizeJ;
        index.m_j = linearIndex - index.m_i*m_sizeJ;

        return index;
    }

    inline std::vector<Index2d> getNeighborhood(Index2d index, int radius = 1, bool vonNeumann = false)
    {
        return getNeighborhood(index.m_i, index.m_j, radius, vonNeumann);
    }

    inline std::vector<Index2d> getNeighborhood(int i, int j, int radius = 1, bool vonNeumann = false)
    {
        std::vector<Index2d> output;
        output.reserve((radius*2+1)*(radius*2+1));
        for(int iOffset = -radius; iOffset <= radius; iOffset++)
        {
            for(int jOffset = -radius; jOffset <= radius; jOffset++)
            {
                if(iOffset == 0 && jOffset == 0) continue;
                if(std::abs(iOffset) + std::abs(jOffset) > radius && vonNeumann) continue;
                Index2d index = Index2d(i + iOffset,j + jOffset);
                if(index.m_i >= 0 && index.m_i < m_sizeI && index.m_j >= 0 && index.m_j < m_sizeJ)
                {
                    output.push_back(index);
                }
            }
        }

        return output;
    }

    inline bool inBounds(Index2d index) const
    {
        return (index.m_i >= 0 && index.m_i < m_sizeI && index.m_j >= 0 && index.m_j < m_sizeJ);
    }

    inline bool inBounds(int i, int j) const
    {
        return (i >= 0 && i < m_sizeI && j >= 0 && j < m_sizeJ);
    }

protected:
    int m_sizeI;
    int m_sizeJ;
};

#endif // LINEARINDEXABLE2D_H
