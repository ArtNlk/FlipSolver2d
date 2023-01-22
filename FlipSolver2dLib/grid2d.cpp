#include "grid2d.h"
#include "mathfuncs.h"

template<class T>
void Grid2d<T>::swap(Grid2d<T> other)
{
    std::swap(m_sizeI,other.m_sizeI);
    std::swap(m_sizeJ,other.m_sizeJ);
    m_data.swap(other.m_data);
}

template<class T>
OOBStrategy &Grid2d<T>::oobStrat()
{
    return m_oobStrat;
}

template<class T>
template<class U,typename std::enable_if<std::is_floating_point<U>::value>::type*>
T Grid2d<T>::interpolateAt(Vertex position)
{
    return simmath::lerpCenteredGrid(position,*this);
}

template<class T>
template<class U,typename std::enable_if<std::is_floating_point<U>::value>::type*>
T Grid2d<T>::interpolateAt(float i, float j)
{
    return simmath::lerpCenteredGrid(i,j,*this);
}

template<class T>
typename std::vector<T>::reference Grid2d<T>::at(int i, int j)
{
    ASSERT_BETWEEN(i,-1,m_sizeI);
    ASSERT_BETWEEN(j,-1,m_sizeJ);
    return m_data[linearIndex(i,j)];
}

template<class T>
typename std::vector<T>::reference Grid2d<T>::at(Index2d index)
{
    ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
    ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
    return m_data[linearIndex(index)];
}

template<class T>
const typename std::vector<T>::const_reference Grid2d<T>::at(int i, int j) const
{
    ASSERT_BETWEEN(i,-1,m_sizeI);
    ASSERT_BETWEEN(j,-1,m_sizeJ);
    return m_data[linearIndex(i,j)];
}

template<class T>
const typename std::vector<T>::const_reference Grid2d<T>::at(Index2d index) const
{
    ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
    ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
    return m_data[linearIndex(index)];
}

template<class T>
void Grid2d<T>::setAt(int i, int j, T value)
{
    ASSERT_BETWEEN(i,-1,m_sizeI);
    ASSERT_BETWEEN(j,-1,m_sizeJ);
    m_data[linearIndex(i,j)] = value;
}

template<class T>
T Grid2d<T>::getAt(Index2d index) const
{
    return getAt(index.m_i, index.m_j);
}

template<class T>
T Grid2d<T>::getAt(int i, int j) const
{
    switch(m_oobStrat)
    {
    case OOB_EXTEND:
        i = std::clamp(i, 0,m_sizeI-1);
        j = std::clamp(j,0,m_sizeJ -1);
        break;
    case OOB_CONST:
        if(!inBounds(i,j))
        {
            return m_oobConst;
        }
        break;
    case OOB_ERROR:
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        break;
    }
    return m_data[linearIndex(i,j)];
}

template<class T>
void Grid2d<T>::fill(T value)
{
    for(int i = 0; i < m_data.size(); i++)
    {
        m_data[i] = value;
    }
}

template<class T>
std::vector<T> &Grid2d<T>::data()
{
    return m_data;
}

template<class T>
void Grid2d<T>::fillRect(T value, Index2d topLeft, Index2d bottomRight)
{
    fillRect(value,topLeft.m_i,topLeft.m_j,bottomRight.m_i,bottomRight.m_j);
}

template<class T>
void Grid2d<T>::fillRect(T value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
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
