#ifndef GRID2D_H
#define GRID2D_H

#include <vector>

#include "linearindexable2d.h"

template<class T>
class Grid2d : public LinearIndexable2d
{
public:
    Grid2d(int sizeI, int sizeJ, T initValue = T()) :
        LinearIndexable2d(sizeI,sizeJ)
    {
        m_data.assign(sizeI*sizeJ,initValue);
    }

    inline T& at(int i, int j)
    {
        ASSERT_BETWEEN(i,0,m_sizeI);
        ASSERT_BETWEEN(j,0,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    inline T& at(Index2d index)
    {
        ASSERT_BETWEEN(i,0,m_sizeI);
        ASSERT_BETWEEN(j,0,m_sizeJ);
        return m_data[linearIndex(index)];
    }

    inline void setAt(int i, int j, T value)
    {
        ASSERT_BETWEEN(i,0,m_sizeI);
        ASSERT_BETWEEN(j,0,m_sizeJ);
        m_data[linearIndex(i,j)] = value;
    }

    inline T getAt(int i, int j) const
    {
        ASSERT_BETWEEN(i,0,m_sizeI);
        ASSERT_BETWEEN(j,0,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    inline void fill(T value)
    {
        for(int i = 0; i < m_data.size(); i++)
        {
            m_data[i] = value;
        }
    }

    inline std::vector<T> &data()
    {
        return m_data;
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

protected:
    std::vector<T> m_data;
};

#endif // GRID2D_H
