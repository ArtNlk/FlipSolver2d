#ifndef MATERIALGRID_H
#define MATERIALGRID_H

#include "grid2d.h"

enum FluidMaterial : int8_t
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
    MaterialGrid(ssize_t sizeI, ssize_t sizeJ, FluidMaterial oobMaterial);

    bool isFluid(Index2d idx) const;

    bool isStrictFluid(Index2d idx) const;

    bool isSolid(Index2d idx) const;

    bool isEmpty(Index2d idx) const;

    bool isSource(Index2d idx) const;

    bool isSink(Index2d idx) const;

    bool isFluid(ssize_t i, ssize_t j) const;

    bool isStrictFluid(ssize_t i, ssize_t j) const;

    bool isSolid(ssize_t i, ssize_t j) const;

    bool isEmpty(ssize_t i, ssize_t j) const;

    bool isSource(ssize_t i, ssize_t j) const;

    bool isSink(ssize_t i, ssize_t j) const;

    bool isFluid(ssize_t i) const;

    bool isStrictFluid(ssize_t i) const;

    bool isSolid(ssize_t i) const;

    bool isEmpty(ssize_t i) const;

    bool isSource(ssize_t i) const;

    bool isSink(ssize_t i) const;

    bool isFluid(size_t i) const;

    bool isStrictFluid(size_t i) const;

    bool isSolid(size_t i) const;

    bool isEmpty(size_t i) const;

    bool isSource(size_t i) const;

    bool isSink(size_t i) const;

    int nonsolidNeighborCount(ssize_t linIdx);

    int nonsolidNeighborCount(ssize_t i, ssize_t j);

    bool uVelocitySampleInside(ssize_t i, ssize_t j) const;

    bool vVelocitySampleInside(ssize_t i, ssize_t j) const;

    bool uSampleAffectedBySolid(ssize_t i, ssize_t j) const;

    bool vSampleAffectedBySolid(ssize_t i, ssize_t j) const;

    VelocitySampleState uVelocitySampleState(ssize_t i, ssize_t j) const;

    VelocitySampleState vVelocitySampleState(ssize_t i, ssize_t j) const;
};

#endif // MATERIALGRID_H
