#ifndef FLUIDGRID_H
#define FLUIDGRID_H

#include <unordered_map>

#include "fluidcell.h"
#include "grid2d.h"
#include "customassert.h"
#include "index2d.h"
#include "logger.h"

class MACFluidGrid : public LinearIndexable2d
{
public:
    MACFluidGrid(int sizeI, int sizeJ);

    void fill(FluidCellMaterial m = FluidCellMaterial::EMPTY,
                     float velocityU = 0.f,
                     float velocityV = 0.f,
                     bool knownFlagU = false,
                     bool knownFlagV = false);

    void fillMaterial(FluidCellMaterial m);

    void fillMaterialRect(FluidCellMaterial value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    void fillMaterialRect(FluidCellMaterial value, Index2d topLeft, Index2d bottomRight);

    void fillVelocityU(float value);

    void fillVelocityURect(double value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    void fillVelocityURect(double value, Index2d topLeft, Index2d bottomRight);

    void fillVelocityV(float value);

    void fillVelocityVRect(double value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    void fillVelocityVRect(double value, Index2d topLeft, Index2d bottomRight);

    void fillKnownFlagsU(bool value);

    void fillKnownFlagURect(bool value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    void fillKnownFlagURect(bool value, Index2d topLeft, Index2d bottomRight);

    void fillKnownFlagsV(bool value);

    void fillKnownFlagVRect(bool value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    void fillKnownFlagVRect(bool value, Index2d topLeft, Index2d bottomRight);

    FluidCellMaterial getMaterial(Index2d index) const;

    FluidCellMaterial getMaterial(int i, int j) const;

    void setMaterial(Index2d index, FluidCellMaterial m);

    void setMaterial(int i, int j, FluidCellMaterial m);

    void setU(Index2d index, double value, bool knownStatus);

    void setU(Index2d index, double value);

    void setU(int i, int j, double value, bool knownStatus);

    void setU(int i, int j, double value);

    void setV(Index2d index, double value, bool knownStatus);

    void setV(Index2d index, double value);

    void setV(int i, int j, double value, bool knownStatus);

    void setV(int i, int j, double value);

    double getU(Index2d index) const;

    double getU(int i, int j) const;

    double getV(Index2d index) const;

    double getV(int i, int j) const;

    void getSize(int& sizeI, int& sizeJ) const;

    int cellCount() const;

    int fluidCellCount() const;

    float trueU(int i, int j);

    float trueV(int i, int j);

//    inline FluidCell at(int i, int j)
//    {
//        return FluidCell(m_materialGrid.at(i,j),
//                         m_velocityGridU.at(i,j),
//                         m_velocityGridV.at(i,j),
//                         m_knownFlagsU.at(i,j),
//                         m_knownFlagsV.at(i,j));
//    }

//    inline FluidCell at(Index2d index)
//    {
//        return FluidCell(m_materialGrid.at(index),
//                         m_velocityGridU.at(index),
//                         m_velocityGridV.at(index),
//                         m_knownFlagsU.at(index),
//                         m_knownFlagsV.at(index));
//    }

    const Grid2d<FluidCellMaterial> &materialGrid();

    Grid2d<float> &velocityGridU();

    Grid2d<float> &velocityGridV();

    Grid2d<bool> &knownFlagsGridU();

    Grid2d<bool> &knownFlagsGridV();

    void updateLinearToFluidMapping();

    int linearFluidIndex(int i, int j);

protected:
    Grid2d<FluidCellMaterial> m_materialGrid;
    Grid2d<float> m_velocityGridU;
    Grid2d<float> m_velocityGridV;
    Grid2d<bool> m_knownFlagsU;
    Grid2d<bool> m_knownFlagsV;
    std::unordered_map<int,int> m_linearToFluidCellIndexMap;
    int m_fluidCellCount;
};

#endif // FLUIDGRID_H
