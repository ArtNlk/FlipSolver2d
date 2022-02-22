#ifndef FLUIDGRID_H
#define FLUIDGRID_H

#include "fluidcell.h"
#include "grid2d.h"
#include "customassert.h"
#include "index2d.h"

class MACFluidGrid : public Grid2d<FluidCell>
{
public:
    MACFluidGrid(int sizeI, int sizeJ) :
        Grid2d(sizeI,sizeJ)
    {
    }

    inline FluidCellMaterial getMaterial(Index2d index) const
    {
        return m_data[linearIndex(index)].getMaterial();
    }

    inline FluidCellMaterial getMaterial(int i, int j) const
    {
        return m_data[linearIndex(i, j)].getMaterial();
    }

    inline void setMaterial(Index2d index, FluidCellMaterial m)
    {
        m_data[linearIndex(index)].setMatrial(m);
    }

    inline void setU(Index2d index, double value)
    {
        m_data[linearIndex(index)].setU(value);
    }

    inline void setV(Index2d index, double value)
    {
        m_data[linearIndex(index)].setV(value);
    }

    inline double getU(Index2d index) const
    {
        return m_data[linearIndex(index)].getU();
    }

    inline double getV(Index2d index) const
    {
        return m_data[linearIndex(index)].getV();
    }

    inline void getSize(int& sizeI, int& sizeJ) const
    {
        sizeI = m_sizeI;
        sizeJ = m_sizeJ;
    }

};

#endif // FLUIDGRID_H
