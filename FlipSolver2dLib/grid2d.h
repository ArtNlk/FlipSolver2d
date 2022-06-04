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

    inline void swap(Grid2d<T> other)
    {
        std::swap(m_sizeI,other.m_sizeI);
        std::swap(m_sizeJ,other.m_sizeJ);
        m_data.swap(other.m_data);
    }

    inline T& at(int i, int j)
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    inline T& at(Index2d index)
    {
        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
        return m_data[linearIndex(index)];
    }

    inline const T& at(int i, int j) const
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    inline const T& at(Index2d index) const
    {
        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
        return m_data[linearIndex(index)];
    }

    inline void setAt(int i, int j, T value)
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        m_data[linearIndex(i,j)] = value;
    }

    inline T getAt(int i, int j) const
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
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

    inline void fillRect(T value, Index2d topLeft, Index2d bottomRight)
    {
        fillRect(value,topLeft.m_i,topLeft.m_j,bottomRight.m_i,bottomRight.m_j);
    }

    inline void fillRect(T value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
    {
        topLeftX = std::max(0,topLeftX);
        topLeftY = std::max(0,topLeftY);
        bottomRightX = std::min(m_sizeI - 1,bottomRightX);
        bottomRightY = std::min(m_sizeJ - 1,bottomRightY);
        for(int i = topLeftX; i <= bottomRightX; i++)
        {
            for(int j = topLeftY; j <= bottomRightY; j++)
            {
                m_data[linearIndex(i,j)] = value;
            }
        }
    }

protected:
    std::vector<T> m_data;
};

template<>
class Grid2d<bool> : public LinearIndexable2d
{
public:
    Grid2d(int sizeI, int sizeJ, bool initValue = false) :
        LinearIndexable2d(sizeI,sizeJ)
    {
        m_data.assign(sizeI*sizeJ,initValue);
    }

    inline std::vector<bool>::reference at(int i, int j)
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    inline std::vector<bool>::reference at(Index2d index)
    {
        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
        return m_data[linearIndex(index)];
    }

    inline std::vector<bool>::const_reference at(int i, int j) const
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    inline std::vector<bool>::const_reference at(Index2d index) const
    {
        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
        return m_data[linearIndex(index)];
    }

    inline void setAt(int i, int j, bool value)
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        m_data[linearIndex(i,j)] = value;
    }

    inline bool getAt(int i, int j) const
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    inline void fill(bool value)
    {
        for(int i = 0; i < m_data.size(); i++)
        {
            m_data[i] = value;
        }
    }

    inline std::vector<bool> &data()
    {
        return m_data;
    }

    inline void fillRect(bool value, Index2d topLeft, Index2d bottomRight)
    {
        fillRect(value,topLeft.m_i,topLeft.m_j,topLeft.m_i,topLeft.m_j);
    }

    inline void fillRect(bool value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
    {
        topLeftX = std::max(0,topLeftX);
        topLeftY = std::max(0,topLeftY);
        bottomRightX = std::min(m_sizeI,bottomRightX);
        bottomRightY = std::min(m_sizeJ,bottomRightY);
        for(int i = topLeftX; i <= bottomRightX; i++)
        {
            for(int j = topLeftY; j <= bottomRightY; j++)
            {
                m_data[linearIndex(i,j)] = value;
            }
        }
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
    std::vector<bool> m_data;
};

#endif // GRID2D_H
