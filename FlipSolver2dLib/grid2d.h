#ifndef GRID2D_H
#define GRID2D_H

#include <algorithm>
#include <vector>

#include "geometry2d.h"
#include "hwinfo.h"
#include "linearindexable2d.h"
#include "mathfuncs.h"
#include "grid_sse42.h"

enum OOBStrategy : char {OOB_EXTEND,
                        OOB_CONST,
                        OOB_ERROR};

template<class T>
class Grid2d : public LinearIndexable2d
{
public:
    Grid2d(int sizeI, int sizeJ,
           T initValue = T(), OOBStrategy oobStrat = OOB_ERROR,
           T oobVal = T(), Vertex gridOffset = Vertex(0.f,0.f)) :
        LinearIndexable2d(sizeI,sizeJ),
        m_cubicInterpFunc(&Grid2d::cubicInterp),
        m_oobStrat(oobStrat),
        m_data(sizeI*sizeJ),
        m_oobConst(oobVal),
        m_gridOffset(gridOffset)
    {
        m_data.assign(sizeI*sizeJ,initValue);
        //initSimdPtrs();
    }

    friend Grid_sse42;

    void swap(Grid2d<T> other)
    {
        std::swap(m_sizeI,other.m_sizeI);
        std::swap(m_sizeJ,other.m_sizeJ);
        m_data.swap(other.m_data);
    }

    OOBStrategy &oobStrat()
    {
        return m_oobStrat;
    }

    template<class U = T, typename std::enable_if<std::is_floating_point<U>::value>::type* = nullptr>
    T interpolateAt(Vertex position) const
    {
        return interpolateAt(position.x(), position.y());
    }

    template<class U = T, typename std::enable_if<std::is_floating_point<U>::value>::type* = nullptr>
    T interpolateAt(float i, float j) const
    {
        //return simmath::lerpCenteredGrid(i, j, *this, m_gridOffset);
        return m_cubicInterpFunc(this,i, j);
    }

    typename std::vector<T>::reference at(int i, int j)
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    typename std::vector<T>::reference at(Index2d index)
    {
        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
        return m_data[linearIndex(index)];
    }

    const typename std::vector<T>::const_reference at(int i, int j) const
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        return m_data[linearIndex(i,j)];
    }

    const typename std::vector<T>::const_reference at(Index2d index) const
    {
        ASSERT_BETWEEN(index.m_i,-1,m_sizeI);
        ASSERT_BETWEEN(index.m_j,-1,m_sizeJ);
        return m_data[linearIndex(index)];
    }

    void setAt(Index2d idx, T value)
    {
        setAt(idx, value);
    }

    void setAt(int i, int j, T value)
    {
        ASSERT_BETWEEN(i,-1,m_sizeI);
        ASSERT_BETWEEN(j,-1,m_sizeJ);
        m_data[linearIndex(i,j)] = value;
    }

    T getAt(Index2d index) const
    {
        return getAt(index.m_i, index.m_j);
    }

    T getAt(int i, int j) const
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
        int idx = linearIndex(i,j);
        return m_data[idx];
    }

    void fill(T value)
    {
        m_data.assign(m_data.size(),value);
    }

    T oobVal() const
    {
        return m_oobConst;
    }

    std::vector<T> &data()
    {
        return m_data;
    }

    const std::vector<T> &data() const
    {
        return m_data;
    }

    void fillRect(T value, Index2d topLeft, Index2d bottomRight)
    {
        fillRect(value,topLeft.m_i,topLeft.m_j,bottomRight.m_i,bottomRight.m_j);
    }

    void fillRect(T value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
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
    template<class U = T, typename std::enable_if<std::is_floating_point<U>::value>::type* = nullptr>
    T static cubicInterp(const Grid2d<T>* grid,float i, float j)
    {
        i+= grid->m_gridOffset.x();
        j+= grid->m_gridOffset.y();

        i -= 0.5;
        j -= 0.5;
        float iFactor = simmath::frac(i);
        float jFactor = simmath::frac(j);
        int iInt = simmath::integr(i);
        int jInt = simmath::integr(j);

        auto calcWeights = [](float f)
        {
            float fSqrd = f*f;
            float fCubed = fSqrd*f;
            std::array<float,4> out = {0.f};
            out[0] = -(1.f/3.f)*f + (1.f/2.f)*fSqrd - (1.f/6.f)*fCubed;
            out[1] = 1.f - fSqrd + (1.f/2.f)*(fCubed-f);
            out[2] = f + (1.f/2.f)*(fSqrd - fCubed);
            out[3] = (1.f/6.f)*(fCubed - f);
            return out;
        };

        std::array<float,4> iWeights = calcWeights(iFactor);
        std::array<float,4> jWeights = calcWeights(jFactor);
        std::array<float,4> temp = {0.f};

        for(int index = 0; index < 4; index++)
        {
            int jCoord = jInt + index - 1;
            temp[index] = iWeights[0]*grid->getAt(i-1,jCoord)
                          + iWeights[1]*grid->getAt(i,jCoord)
                          + iWeights[2]*grid->getAt(i+1,jCoord)
                          + iWeights[3]*grid->getAt(i+2,jCoord);
        }

        float output = jWeights[0]*temp[0]
                       + jWeights[1]*temp[1]
                       + jWeights[2]*temp[2]
                       + jWeights[3]*temp[3];
        return output;
    }

    //Simds for float
    template<class U = T, typename std::enable_if<std::is_same<U,float>::value>::type* = nullptr>
    void initSimdPtrs()
    {
        switch(HwInfo::i().getSimdLevel())
        {
        case SIMD_LEVEL_NONE:
            break;
        case SIMD_LEVEL_SSE42:
            m_cubicInterpFunc = &Grid_sse42::cubicInterpF;
            break;
        case SIMD_LEVEL_SSE4a_XOP_FMA:
        case SIMD_LEVEL_AVX:
        case SIMD_LEVEL_AVX2:
        case SIMD_LEVEL_AVX512:
            break;
        }
    }

    //Simds for double
    template<class U = T, typename std::enable_if<std::is_same<U,double>::value>::type* = nullptr>
    void initSimdPtrs()
    {
        switch(HwInfo::i().getSimdLevel())
        {
        case SIMD_LEVEL_NONE:
            break;
        case SIMD_LEVEL_SSE42:
            //m_cubicInterpFunc = &Grid_sse42::cubicInterpD;
            break;
        case SIMD_LEVEL_SSE4a_XOP_FMA:
        case SIMD_LEVEL_AVX:
        case SIMD_LEVEL_AVX2:
        case SIMD_LEVEL_AVX512:
            break;
        }
    }

    //Cubic interp placeholder for non-float types
    template<class U = T, typename std::enable_if<!std::is_floating_point<U>::value>::type* = nullptr>
    T static cubicInterp(const Grid2d<T>* grid,float i, float j)
    {
        return grid->getAt(i,j);
    }

    template<class U = T, typename std::enable_if<!std::is_floating_point<U>::value>::type* = nullptr>
    void initSimdPtrs()
    {
        return;
    }

protected:
    std::vector<T> m_data;
    OOBStrategy m_oobStrat;
    Vertex m_gridOffset;
    T m_oobConst;

    T(*m_cubicInterpFunc)(const Grid2d*,float,float);
};

#endif // GRID2D_H
