#ifndef FLUIDGRID_H
#define FLUIDGRID_H

#include <unordered_map>

#include "grid2d.h"
#include "materialgrid.h"
#include "customassert.h"
#include "index2d.h"
#include "logger.h"
#include "geometry2d.h"
#include "staggeredvelocitygrid.h"

class MACFluidGrid : public LinearIndexable2d
{
public:
    MACFluidGrid(int sizeI, int sizeJ);

    void fill(FluidMaterial m = FluidMaterial::EMPTY,
                     float velocityU = 0.f,
                     float velocityV = 0.f,
                     bool knownFlagU = false,
                     bool knownFlagV = false);

    void fillMaterial(FluidMaterial m);

    void fillMaterialRect(FluidMaterial value, int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    void fillMaterialRect(FluidMaterial value, Index2d topLeft, Index2d bottomRight);

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

    FluidMaterial getMaterial(Index2d index) const;

    FluidMaterial getMaterial(int i, int j) const;

    void setMaterial(Index2d index, FluidMaterial m);

    void setMaterial(int i, int j, FluidMaterial m);



    void setFluidU(Index2d index, float value, bool knownStatus);

    void setFluidU(Index2d index, float value);

    void setFluidU(int i, int j, float value, bool knownStatus);

    void setFluidU(int i, int j, float value);

    void setFluidV(Index2d index, float value, bool knownStatus);

    void setFluidV(Index2d index, float value);

    void setFluidV(int i, int j, float value, bool knownStatus);

    void setFluidV(int i, int j, float value);

    void setAirU(int i, int j, float value);

    void setAirV(int i, int j, float value);

    void getSize(int& sizeI, int& sizeJ) const;

    float getFluidU(int i, int j);

    float getFluidV(int i, int j);

    float getAirU(int i, int j);

    float getAirV(int i, int j);

    int cellCount() const;

    int fluidCellCount() const;

    Vertex fluidVelocityAt(float i, float j);

    Vertex fluidVelocityAt(Vertex position);

    Vertex airVelocityAt(float i, float j);

    Vertex airVelocityAt(Vertex position);

    float viscosityAt(Vertex position);

    Grid2d<FluidMaterial> &materialGrid();

    StaggeredVelocityGrid& fluidVelocityGrid();

    StaggeredVelocityGrid& airVelocityGrid();

    StaggeredVelocityGrid& savedFluidVelocityGrid();

    StaggeredVelocityGrid& savedAirVelocityGrid();

    Grid2d<float> &fluidVelocityGridU();

    Grid2d<float> &fluidVelocityGridV();

    Grid2d<float> &airVelocityGridU();

    Grid2d<float> &airVelocityGridV();

    Grid2d<bool> &knownFluidFlagsGridU();

    Grid2d<bool> &knownFluidFlagsGridV();

    Grid2d<bool> &knownAirFlagsGridU();

    Grid2d<bool> &knownAirFlagsGridV();

    Grid2d<bool> &knownFlagsCenteredParams();

    Grid2d<float> &viscosityGrid();

    Grid2d<float> &solidSdfGrid();

    Grid2d<float> &fluidSdfGrid();

    Grid2d<float> &airSdfGrid();

    Grid2d<float> &testGrid();

    Grid2d<float> &temperatureGrid();

    Grid2d<float> &smokeConcentrationGrid();

    Grid2d<float> &divergenceControlGrid();

    Grid2d<int> &fluidParticleCountGrid();

    Grid2d<int> &airParticleCountGrid();

    float solidSdf(int i, int j);

    float solidSdfAt(float i, float j);

    float solidSdf(Index2d index);

    void setSolidSdf(int i, int j, float value);

    void setSolidSdf(Index2d index, float value);

    float fluidSdfAt(float i, float j);

    float viscosity(int i, int j);

    float viscosityAt(float i, float j);

    float viscosity(Index2d index);

    void setViscosity(int i, int j, float value);

    void setViscosity(Index2d index, float value);

    float temperature(int i, int j);

    float temperatureAt(float i, float j);

    float temperatureAt(Vertex& v);

    void setTemperature(int i, int j, float value);

    float smokeConcentration(int i, int j);

    float smokeConcentrationAt(float i, float j);

    float smokeConcentrationAt(Vertex v);

    void setSmokeConcentration(int i, int j, float value);

    float &divergenceControl(int i, int j);

    void setDivergeceControl(int i, int j, float value);

    int closestSolidId(int i, int j);

    std::vector<int> nearValidSolidIds(int i, int j);

    Vertex closestSolidSurfacePoint(Vertex pos);

    Vertex closesFluidSurfacePoint(Vertex pos);

    void setEmitterId(int i, int j, int id);

    int emitterId(int i, int j);

    void setSolidId(int i, int j, int solidId);

    int solidId(int i, int j);

    float getFaceFractionUSample(int i, int j);

    float getFaceFractionVSample(int i, int j);

protected:  


    StaggeredVelocityGrid m_airVelocityGrid;
    StaggeredVelocityGrid m_savedAirVelocityGrid;
    Grid2d<float> m_temperature;
    Grid2d<float> m_smokeConcentration;
    Grid2d<int> m_airParticleCounts;
    Grid2d<float> m_airSdf;
    int m_validUVelocitySampleCount;
    int m_validVVelocitySampleCount;
    int m_fluidCellCount;
};

#endif // FLUIDGRID_H
