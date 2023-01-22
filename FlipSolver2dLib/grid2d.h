#ifndef GRID2D_H
#define GRID2D_H

#include <algorithm>
#include <vector>

#include "geometry2d.h"
#include "linearindexable2d.h"

enum OOBStrategy : char {OOB_EXTEND,
                        OOB_CONST,
                        OOB_ERROR};

template<class T>
class Grid2d : public LinearIndexable2d
{
public:
    Grid2d(int sizeI, int sizeJ, T initValue = T(), OOBStrategy oobStrat = OOB_ERROR, T oobVal = T()) :
        LinearIndexable2d(sizeI,sizeJ),
        m_oobStrat(oobStrat),
        m_oobConst(oobVal)
    {
        m_data.assign(sizeI*sizeJ,initValue);
    }

    void swap(Grid2d<T> other);

    OOBStrategy &oobStrat();

    template<class U = T, typename std::enable_if<std::is_floating_point<U>::value>::type* = nullptr>
    T interpolateAt(Vertex position);

    template<class U = T, typename std::enable_if<std::is_floating_point<U>::value>::type* = nullptr>
    T interpolateAt(float i, float j);

    typename std::vector<T>::reference at(int i, int j);

    typename std::vector<T>::reference at(Index2d index);

    const typename std::vector<T>::const_reference at(int i, int j) const;

    const typename std::vector<T>::const_reference at(Index2d index) const;

    void setAt(int i, int j, T value);

    T getAt(Index2d index) const;

    T getAt(int i, int j) const;

    void fill(T value);

    std::vector<T> &data();

    void fillRect(T value, Index2d topLeft, Index2d bottomRight);

    void fillRect(T value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

protected:
    std::vector<T> m_data;
    OOBStrategy m_oobStrat;
    T m_oobConst;
};

//template<>
//class Grid2d<bool> : public LinearIndexable2d
//{
//public:
//    Grid2d(int sizeI, int sizeJ, bool initValue = false, OOBStrategy oobStrat = OOB_ERROR, bool oobVal = false) :
//        LinearIndexable2d(sizeI,sizeJ),
//        m_oobStrat(oobStrat),
//        m_oobConst(oobVal)
//    {
//        m_data.assign(sizeI*sizeJ,initValue);
//    }

//    inline OOBStrategy &oobStrat()
//    {
//        return m_oobStrat;
//    }

//    inline std::vector<bool>::reference at(int i, int j)
//    {
//        ASSERT_BETWEEN(i,-1,m_sizeI);
//        ASSERT_BETWEEN(j,-1,m_sizeJ);
//        return m_data[linearIndex(i,j)];
//    }

//    inline std::vector<bool>::reference at(Index2d index)
//    {
//        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
//        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
//        return m_data[linearIndex(index)];
//    }

//    inline std::vector<bool>::const_reference at(int i, int j) const
//    {
//        ASSERT_BETWEEN(i,-1,m_sizeI);
//        ASSERT_BETWEEN(j,-1,m_sizeJ);
//        return m_data[linearIndex(i,j)];
//    }

//    inline std::vector<bool>::const_reference at(Index2d index) const
//    {
//        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
//        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
//        return m_data[linearIndex(index)];
//    }

//    inline void setAt(int i, int j, bool value)
//    {
//        ASSERT_BETWEEN(i,-1,m_sizeI);
//        ASSERT_BETWEEN(j,-1,m_sizeJ);
//        m_data[linearIndex(i,j)] = value;
//    }

//    inline bool getAt(Index2d index) const
//    {
//        return getAt(index.m_i, index.m_j);
//    }

//    inline bool getAt(int i, int j) const
//    {
//        switch(m_oobStrat)
//        {
//        case OOB_EXTEND:
//            i = std::clamp(i, 0,m_sizeI-1);
//            j = std::clamp(j,0,m_sizeJ -1);
//            break;
//        case OOB_CONST:
//            if(!inBounds(i,j))
//            {
//                return m_oobConst;
//            }
//            break;
//        case OOB_ERROR:
//            ASSERT_BETWEEN(i,-1,m_sizeI);
//            ASSERT_BETWEEN(j,-1,m_sizeJ);
//            break;
//        }
//        return m_data[linearIndex(i,j)];
//    }

//    inline void fill(bool value)
//    {
//        for(int i = 0; i < m_data.size(); i++)
//        {
//            m_data[i] = value;
//        }
//    }

//    inline std::vector<bool> &data()
//    {
//        return m_data;
//    }

//    inline void fillRect(bool value, Index2d topLeft, Index2d bottomRight)
//    {
//        fillRect(value,topLeft.m_i,topLeft.m_j,topLeft.m_i,topLeft.m_j);
//    }

//    inline void fillRect(bool value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
//    {
//        topLeftX = std::max(0,topLeftX);
//        topLeftY = std::max(0,topLeftY);
//        bottomRightX = std::min(m_sizeI,bottomRightX);
//        bottomRightY = std::min(m_sizeJ,bottomRightY);
//        for(int i = topLeftX; i <= bottomRightX; i++)
//        {
//            for(int j = topLeftY; j <= bottomRightY; j++)
//            {
//                m_data[linearIndex(i,j)] = value;
//            }
//        }
//    }

//    inline std::vector<Index2d> getNeighborhood(Index2d index, int radius = 1, bool vonNeumann = false)
//    {
//        return getNeighborhood(index.m_i, index.m_j, radius, vonNeumann);
//    }

//    inline std::vector<Index2d> getNeighborhood(int i, int j, int radius = 1, bool vonNeumann = false)
//    {
//        std::vector<Index2d> output;
//        output.reserve((radius*2+1)*(radius*2+1));
//        for(int iOffset = -radius; iOffset <= radius; iOffset++)
//        {
//            for(int jOffset = -radius; jOffset <= radius; jOffset++)
//            {
//                if(iOffset == 0 && jOffset == 0) continue;
//                if(std::abs(iOffset) + std::abs(jOffset) > radius && vonNeumann) continue;
//                Index2d index = Index2d(i + iOffset,j + jOffset);
//                if(index.m_i >= 0 && index.m_i < m_sizeI && index.m_j >= 0 && index.m_j < m_sizeJ)
//                {
//                    output.push_back(index);
//                }
//            }
//        }

//        return output;
//    }

//protected:
//    std::vector<bool> m_data;
//    OOBStrategy m_oobStrat;
//    bool m_oobConst;
//};

#endif // GRID2D_H
