#ifndef FLUIDGRID_H
#define FLUIDGRID_H

#include <unordered_map>

#include "fluidcell.h"
#include "grid2d.h"
#include "customassert.h"
#include "index2d.h"
#include "logger.h"
#include "geometry2d.h"

struct PairHash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    std::size_t lhs = std::hash<T>()(x.first);
    std::size_t rhs = std::hash<U>()(x.second);
    return rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  }
};

//Bit mask:
//1<<0 - partial submersion?
//1<<1 - Fluid
//1<<2 - Solid
//1<<3 - Air
enum VelocitySampleState : char {IN_FLUID = 0b0100,
                                 IN_SOLID = 0b0010,
                                 IN_AIR = 0b0001,
                                 BORDER_FLUID_SOLID = 0b1110,
                                 BORDER_FLUID_AIR = 0b1101,
                                 BORDER_SOLID_AIR = 0b1011};

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

    bool uVelocitySampleInside(int i, int j);

    bool vVelocitySampleInside(int i, int j);

    bool uSampleAffectedBySolid(int i, int j);

    bool vSampleAffectedBySolid(int i, int j);

    VelocitySampleState uVelocitySampleState(int i, int j);

    VelocitySampleState vVelocitySampleState(int i, int j);

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

    float viscosityAt(Vertex position);

    Grid2d<FluidCellMaterial> &materialGrid();

    Grid2d<float> &velocityGridU();

    Grid2d<float> &velocityGridV();

    Grid2d<bool> &knownFlagsGridU();

    Grid2d<bool> &knownFlagsGridV();

    Grid2d<bool> &knownFlagsCenteredParams();

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

    int closestSolidId(int i, int j);

    std::vector<int> nearValidSolidIds(int i, int j);

    Vertex closestSurfacePoint(Vertex pos);

    void updateLinearFluidViscosityMapping();

    int linearViscosityVelocitySampleIndexU(int i, int j);

    int linearViscosityVelocitySampleIndexV(int i, int j);

    int validVelocitySampleCountU();

    int validVelocitySampleCountV();

    void setEmitterId(int i, int j, int id);

    int emitterId(int i, int j);

    void setSolidId(int i, int j, int solidId);

    int solidId(int i, int j);

protected:  
    void updateValidULinearMapping();

    void updateValidVLinearMapping();

    Grid2d<FluidCellMaterial> m_materialGrid;
    Grid2d<float> m_velocityGridU;
    Grid2d<float> m_velocityGridV;
    Grid2d<float> m_sdf;
    Grid2d<bool> m_knownFlagsU;
    Grid2d<bool> m_knownFlagsV;
    Grid2d<bool> m_knownCenteredParams;
    Grid2d<float> m_viscosityGrid;
    Grid2d<int> m_emitterId;
    Grid2d<int> m_solidId;
    std::unordered_map<std::pair<int,int>,int,PairHash> m_uVelocitySamplesMap;
    std::unordered_map<std::pair<int,int>,int,PairHash> m_vVelocitySamplesMap;
    int m_validUVelocitySampleCount;
    int m_validVVelocitySampleCount;
    int m_fluidCellCount;
};

#endif // FLUIDGRID_H
