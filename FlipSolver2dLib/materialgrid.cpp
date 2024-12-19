#include "materialgrid.h"
#include "grid2d.h"
#include "index2d.h"

MaterialGrid::MaterialGrid(ssize_t sizeI, ssize_t sizeJ, FluidMaterial oobMaterial) :
    Grid2d(sizeI, sizeJ, FluidMaterial::EMPTY, OOB_EXTEND, oobMaterial)
{
}

bool MaterialGrid::isFluid(Index2d idx) const
{
    return isFluid(idx.i, idx.j);
}

bool MaterialGrid::isStrictFluid(Index2d idx) const
{
    return isStrictFluid(idx.i, idx.j);
}

bool MaterialGrid::isSolid(Index2d idx) const
{
    return isSolid(idx.i, idx.j);
}

bool MaterialGrid::isEmpty(Index2d idx) const
{
    return isEmpty(idx.i, idx.j);
}

bool MaterialGrid::isSource(Index2d idx) const
{
    return isSource(idx.i, idx.j);
}

bool MaterialGrid::isSink(Index2d idx) const
{
    return isSink(idx.i, idx.j);
}

bool MaterialGrid::isFluid(ssize_t i, ssize_t j) const
{
    return fluidTest(getAt(i,j));
}

bool MaterialGrid::isStrictFluid(ssize_t i, ssize_t j) const
{
    return strictFluidTest(getAt(i,j));
}

bool MaterialGrid::isSolid(ssize_t i, ssize_t j) const
{
    return solidTest(getAt(i,j));
}

bool MaterialGrid::isEmpty(ssize_t i, ssize_t j) const
{
    return emptyTest(getAt(i,j));
}

bool MaterialGrid::isSource(ssize_t i, ssize_t j) const
{
    return sourceTest(getAt(i,j));
}

bool MaterialGrid::isSink(ssize_t i, ssize_t j) const
{
    return sinkTest(getAt(i,j));
}

bool MaterialGrid::isFluid(ssize_t i) const
{
    return fluidTest(getAt(i));
}

bool MaterialGrid::isStrictFluid(ssize_t i) const
{
    return strictFluidTest(getAt(i));
}

bool MaterialGrid::isSolid(ssize_t i) const
{
    return solidTest(getAt(i));
}

bool MaterialGrid::isEmpty(ssize_t i) const
{
    return emptyTest(getAt(i));
}

bool MaterialGrid::isSource(ssize_t i) const
{
    return sourceTest(getAt(i));
}

bool MaterialGrid::isSink(ssize_t i) const
{
    return sinkTest(getAt(i));
}

bool MaterialGrid::isFluid(size_t i) const
{
    return fluidTest(m_data.at(i));
}

bool MaterialGrid::isStrictFluid(size_t i) const
{
    return strictFluidTest(m_data.at(i));
}

bool MaterialGrid::isSolid(size_t i) const
{
    return solidTest(m_data.at(i));
}

bool MaterialGrid::isEmpty(size_t i) const
{
    return emptyTest(m_data.at(i));
}

bool MaterialGrid::isSource(size_t i) const
{
    return sourceTest(m_data.at(i));
}

bool MaterialGrid::isSink(size_t i) const
{
    return sinkTest(m_data.at(i));
}

int MaterialGrid::nonsolidNeighborCount(ssize_t linIdx)
{
    const Index2d i2d = index2d(linIdx);

    return nonsolidNeighborCount(i2d.i, i2d.j);
}

int MaterialGrid::nonsolidNeighborCount(ssize_t i, ssize_t j)
{
    return !isSolid(i-1,j) + !isSolid(i+1,j) + !isSolid(i,j-1) + !isSolid(i,j+1);
}

bool MaterialGrid::uVelocitySampleInside(ssize_t i, ssize_t j) const
{
    return (!isEmpty(i,j) && !isEmpty(i-1, j));
}

bool MaterialGrid::vVelocitySampleInside(ssize_t i, ssize_t j) const
{
    return (!isEmpty(i,j) && !isEmpty(i, j-1));
}

bool MaterialGrid::uSampleAffectedBySolid(ssize_t i, ssize_t j) const
{
    return (isSolid(i,j+1) || isSolid(i-1, j+1) ||
            isSolid(i,j) || isSolid(i-1, j) ||
            isSolid(i,j-1) || isSolid(i-1, j-1));
}

bool MaterialGrid::vSampleAffectedBySolid(ssize_t i, ssize_t j) const
{
    return (isSolid(i+1,j) || isSolid(i+1, j-1) ||
            isSolid(i,j) || isSolid(i, j-1) ||
            isSolid(i-1,j) || isSolid(i-1, j-1));
}

VelocitySampleState MaterialGrid::uVelocitySampleState(ssize_t i, ssize_t j) const
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

VelocitySampleState MaterialGrid::vVelocitySampleState(ssize_t i, ssize_t j) const
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
