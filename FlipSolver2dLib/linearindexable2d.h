#ifndef LINEARINDEXABLE2D_H
#define LINEARINDEXABLE2D_H

#include "index2d.h"
#include "customassert.h"

class LinearIndexable2d
{
public:
    LinearIndexable2d(int sizeI, int sizeJ);

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
        ASSERT_BETWEEN(index.i,0,m_sizeJ);
        ASSERT_BETWEEN(index.j,0,m_sizeJ);
        return index.i * m_sizeJ + index.j;
    }

    inline Index2d index2d(int linearIndex) const
    {
        ASSERT_BETWEEN(linearIndex,0,m_sizeJ*m_sizeI);
        Index2d index;
        index.i = linearIndex / m_sizeJ;
        index.j = linearIndex - index.i*m_sizeJ;

        return index;
    }

protected:
    int m_sizeI;
    int m_sizeJ;
};

#endif // LINEARINDEXABLE2D_H
