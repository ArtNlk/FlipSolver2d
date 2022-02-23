#ifndef FLUIDCELL_H
#define FLUIDCELL_H

enum FluidCellMaterial : char
{
    FLUID = 0,
    SOLID = 1,
    EMPTY = 2
};

class FluidCell
{
public:
    FluidCell() :
        m_material(EMPTY),
        m_velocityU(0.0),
        m_velocityV(0.0)
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

    inline void setU(double value)
    {
        m_velocityU = value;
    }

    inline double getU() const
    {
        return m_velocityU;
    }

    inline double& U()
    {
        return m_velocityU;
    }

    inline void setV(double value)
    {
        m_velocityV = value;
    }

    inline double getV() const
    {
        return m_velocityV;
    }

    inline double& V()
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

    inline bool isKnownU() const
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

    inline bool isKnownV() const
    {
        return m_knownV;
    }

protected:
    FluidCellMaterial m_material;
    double m_velocityU;
    double m_velocityV;
    bool m_knownU;
    bool m_knownV;
};

#endif // FLUIDCELL_H
