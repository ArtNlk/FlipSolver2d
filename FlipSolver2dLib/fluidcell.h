#ifndef FLUIDCELL_H
#define FLUIDCELL_H

#include <vector>

enum FluidCellMaterial : char
{
    FLUID = 0b01000000,
    SOURCE =0b01000001,
    SOLID = 0b00100000,
    SINK =  0b00010010,
    EMPTY = 0b00010000
};

#define fluidTest(x) ((x & FluidCellMaterial::FLUID) != 0)
#define emptyTest(x) ((x & FluidCellMaterial::EMPTY) != 0)
#define solidTest(x) ((x & FluidCellMaterial::SOLID) != 0)
#define sourceTest(x) (x == FluidCellMaterial::SOURCE)
#define sinkTest(x) (x == FluidCellMaterial::SINK)

class FluidCell
{
public:
    FluidCell(FluidCellMaterial& material, float& u, float& v, std::vector<bool>::reference uKnown, std::vector<bool>::reference vKnown) :
        m_material(material),
        m_velocityU(u),
        m_velocityV(v),
        m_knownU(uKnown),
        m_knownV(vKnown)
    {
    }

    inline void setMatrial(FluidCellMaterial m)
    {
        m_material = m;
    }

    inline FluidCellMaterial getMaterial() const
    {
        return m_material;
    }

    inline FluidCellMaterial& material()
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
    FluidCellMaterial& m_material;
    float& m_velocityU;
    float& m_velocityV;
    std::vector<bool>::reference m_knownU;
    std::vector<bool>::reference m_knownV;
};

#endif // FLUIDCELL_H
