#include "materialgrid.h"
#include "grid2d.h"

MaterialGrid::MaterialGrid(int sizeI, int sizeJ, FluidMaterial oobMaterial) :
    Grid2d(sizeI, sizeJ, FluidMaterial::EMPTY, OOB_EXTEND, oobMaterial)
{
}

bool MaterialGrid::isFluid(Index2d idx)
{
    return isFluid(idx.m_i, idx.m_j);
}

bool MaterialGrid::isStrictFluid(Index2d idx)
{
    return isStrictFluid(idx.m_i, idx.m_j);
}

bool MaterialGrid::isSolid(Index2d idx)
{
    return isSolid(idx.m_i, idx.m_j);
}

bool MaterialGrid::isEmpty(Index2d idx)
{
    return isEmpty(idx.m_i, idx.m_j);
}

bool MaterialGrid::isSource(Index2d idx)
{
    return isSource(idx.m_i, idx.m_j);
}

bool MaterialGrid::isSink(Index2d idx)
{
    return isSink(idx.m_i, idx.m_j);
}

bool MaterialGrid::isFluid(int i, int j)
{
    return fluidTest(getAt(i,j));
}

bool MaterialGrid::isStrictFluid(int i, int j)
{
    return strictFluidTest(getAt(i,j));
}

bool MaterialGrid::isSolid(int i, int j)
{
    return solidTest(getAt(i,j));
}

bool MaterialGrid::isEmpty(int i, int j)
{
    return emptyTest(getAt(i,j));
}

bool MaterialGrid::isSource(int i, int j)
{
    return sourceTest(getAt(i,j));
}

bool MaterialGrid::isSink(int i, int j)
{
    return sinkTest(getAt(i,j));
}

bool MaterialGrid::uVelocitySampleInside(int i, int j)
{
    return (!isEmpty(i,j) && !isEmpty(i-1, j));
}

bool MaterialGrid::vVelocitySampleInside(int i, int j)
{
    return (!isEmpty(i,j) && !isEmpty(i, j-1));
}

bool MaterialGrid::uSampleAffectedBySolid(int i, int j)
{
    return (isSolid(i,j+1) || isSolid(i-1, j+1) ||
            isSolid(i,j) || isSolid(i-1, j) ||
            isSolid(i,j-1) || isSolid(i-1, j-1));
}

bool MaterialGrid::vSampleAffectedBySolid(int i, int j)
{
    return (isSolid(i+1,j) || isSolid(i+1, j-1) ||
            isSolid(i,j) || isSolid(i, j-1) ||
            isSolid(i-1,j) || isSolid(i-1, j-1));
}

VelocitySampleState MaterialGrid::uVelocitySampleState(int i, int j)
{
    char currentMat = at(i,j) >> 4;
    char neighborMat = at(i-1,j) >> 4;
    char result = 0;

    result |= currentMat;
    result |= neighborMat;
    if ((result & (result - 1)) != 0)
    {
        result &= ~(0b1000);
    }

    return static_cast<VelocitySampleState>(result);
}

VelocitySampleState MaterialGrid::vVelocitySampleState(int i, int j)
{
    char currentMat = at(i,j) >> 4;
    char neighborMat = at(i-1,j) >> 4;
    char result = 0;

    result |= currentMat;
    result |= neighborMat;
    if ((result & (result - 1)) == 0)
    {
        result &= 0b1000;
    }

    return static_cast<VelocitySampleState>(result);
}
