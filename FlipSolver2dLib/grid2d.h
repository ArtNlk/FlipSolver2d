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

protected:
    std::vector<T> m_data;
};

#endif // GRID2D_H
