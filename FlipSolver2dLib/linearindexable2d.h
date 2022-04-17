#ifndef LINEARINDEXABLE2D_H
#define LINEARINDEXABLE2D_H

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
        ASSERT_BETWEEN(i,0,m_sizeJ);
        ASSERT_BETWEEN(j,0,m_sizeJ);
        return i * m_sizeJ + j;
    }

    inline int linearIndex(Index2d index) const
    {
        ASSERT_BETWEEN(index.m_i,0,m_sizeJ);
        ASSERT_BETWEEN(index.m_j,0,m_sizeJ);
        return index.m_i * m_sizeJ + index.m_j;
    }

    inline Index2d index2d(int linearIndex) const
    {
        ASSERT_BETWEEN(linearIndex,0,m_sizeJ*m_sizeI);
        Index2d index;
        index.m_i = linearIndex / m_sizeJ;
        index.m_j = linearIndex - index.m_i*m_sizeJ;

        return index;
    }

protected:
    int m_sizeI;
    int m_sizeJ;
};

#endif // LINEARINDEXABLE2D_H
