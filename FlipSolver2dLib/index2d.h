#ifndef INDEX2D_H
#define INDEX2D_H


struct Index2d
{
    Index2d() :
        m_i(0),
        m_j(0)
    {
    }

    Index2d(int i, int j) :
        m_i(i),
        m_j(j)
    {
    }

    int m_i;
    int m_j;
};

#endif // INDEX2D_H
