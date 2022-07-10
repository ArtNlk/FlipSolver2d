#ifndef FLUIDGRID_H
#define FLUIDGRID_H

#include <unordered_map>

#include "fluidcell.h"
#include "grid2d.h"
#include "customassert.h"
#include "index2d.h"
#include "logger.h"
#include "geometry2d.h"

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

    void fillVelocityURect(float value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    void fillVelocityURect(float value, Index2d topLeft, Index2d bottomRight);

    void fillVelocityV(float value);

    void fillVelocityVRect(float value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    void fillVelocityVRect(float value, Index2d topLeft, Index2d bottomRight);

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

    bool isFluid(int i, int j);

    bool isSolid(int i, int j);

    bool isEmpty(int i, int j);

    bool isSource(int i, int j);

    bool isSink(int i, int j);

    void setU(Index2d index, float value, bool knownStatus);

    void setU(Index2d index, float value);

    void setU(int i, int j, float value, bool knownStatus);

    void setU(int i, int j, float value);

    void setV(Index2d index, float value, bool knownStatus);

    void setV(Index2d index, float value);

    void setV(int i, int j, float value, bool knownStatus);

    void setV(int i, int j, float value);

    float getU(Index2d index) const;

    float getU(int i, int j) const;

    float getV(Index2d index) const;

    float getV(int i, int j) const;

    void getSize(int& sizeI, int& sizeJ) const;

    int cellCount() const;

    int fluidCellCount() const;

    Vertex velocityAt(float i, float j);

    Vertex velocityAt(Vertex position);

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

    Grid2d<FluidCellMaterial> &materialGrid();

    Grid2d<float> &velocityGridU();

    Grid2d<float> &velocityGridV();

    Grid2d<bool> &knownFlagsGridU();

    Grid2d<bool> &knownFlagsGridV();

    Grid2d<float> &viscosityGrid();

    Grid2d<float> &sdfGrid();

    float sdf(int i, int j);

    float sdfAt(float i, float j);

    float sdf(Index2d index);

    void setSdf(int i, int j, float value);

    void setSdf(Index2d index, float value);

    float viscosity(int i, int j);

    float viscosityAt(float i, float j);

    float viscosity(Index2d index);

    void setViscosity(int i, int j, float value);

    void setViscosity(Index2d index, float value);

    Vertex closestSurfacePoint(Vertex pos);

    void updateLinearToFluidMapping();

    int linearFluidIndex(int i, int j);

protected:
    Grid2d<FluidCellMaterial> m_materialGrid;
    Grid2d<float> m_velocityGridU;
    Grid2d<float> m_velocityGridV;
    Grid2d<float> m_sdf;
    Grid2d<bool> m_knownFlagsU;
    Grid2d<bool> m_knownFlagsV;
    Grid2d<float> m_viscosity;
    std::unordered_map<int,int> m_linearToFluidCellIndexMap;
    int m_fluidCellCount;
};

#endif // FLUIDGRID_H
