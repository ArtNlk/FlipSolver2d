#ifndef FLUIDCELL_H
#define FLUIDCELL_H

#include <vector>

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

class FluidCell
{
public:
    FluidCell(FluidMaterial& material, float& u, float& v, std::vector<bool>::reference uKnown, std::vector<bool>::reference vKnown) :
        m_material(material),
        m_velocityU(u),
        m_velocityV(v),
        m_knownU(uKnown),
        m_knownV(vKnown)
    {
    }

    inline void setMatrial(FluidMaterial m)
    {
        m_material = m;
    }

    inline FluidMaterial getMaterial() const
    {
        return m_material;
    }

    inline FluidMaterial& material()
    {
        return m_material;
    }

    inline void setU(float value)
    {
        m_velocityU = value;
    }

    inline float getU() const
    {
        return m_velocityU;
    }

    inline float& U()
    {
        return m_velocityU;
    }

    inline void setV(float value)
    {
        m_velocityV = value;
    }

    inline float getV() const
    {
        return m_velocityV;
    }

    inline float& V()
    {
        return m_velocityV;
    }

    inline void setKnownStatusU(bool status)
    {
        m_knownU = status;
    }

    inline void setKnownU()
    {
        m_knownU = true;
    }

    inline void setUnknownU()
    {
        m_knownU = false;
    }

    inline std::vector<bool>::reference knownU() const
    {
        return m_knownU;
    }

    inline void setKnownStatusV(bool status)
    {
        m_knownV = status;
    }

    inline void setKnownV()
    {
        m_knownV = true;
    }

    inline void setUnknownV()
    {
        m_knownV = false;
    }

    inline std::vector<bool>::reference knownV() const
    {
        return m_knownV;
    }

protected:
    FluidMaterial& m_material;
    float& m_velocityU;
    float& m_velocityV;
    std::vector<bool>::reference m_knownU;
    std::vector<bool>::reference m_knownV;
};

#endif // FLUIDCELL_H
