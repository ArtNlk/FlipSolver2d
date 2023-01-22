#ifndef MATERIALGRID_H
#define MATERIALGRID_H

#include "grid2d.h"

enum FluidMaterial : char
{
    FLUID = 0b01000000,
    SOURCE =0b01000001,
    SOLID = 0b00100000,
    SINK =  0b00010010,
    EMPTY = 0b00010000
};

constexpr bool fluidTest(FluidMaterial m) { return ((m & FluidMaterial::FLUID) != 0);}
constexpr bool strictFluidTest(FluidMaterial m) { return (m == FluidMaterial::FLUID);}
constexpr bool emptyTest(FluidMaterial m) { return ((m & FluidMaterial::EMPTY) != 0);}
constexpr bool solidTest(FluidMaterial m) { return ((m & FluidMaterial::SOLID) != 0);}
constexpr bool sourceTest(FluidMaterial m) { return (m == FluidMaterial::SOURCE);}
constexpr bool sinkTest(FluidMaterial m) { return (m == FluidMaterial::SINK);}

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

class MaterialGrid : public Grid2d<FluidMaterial>
{
public:
    MaterialGrid(int sizeI, int sizeJ, FluidMaterial oobMaterial);

    bool isFluid(int i, int j);

    bool isStrictFluid(int i, int j);

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
};

#endif // MATERIALGRID_H
