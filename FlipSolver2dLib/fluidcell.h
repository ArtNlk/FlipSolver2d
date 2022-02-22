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

protected:
    FluidCellMaterial m_material;
    double m_velocityU;
    double m_velocityV;
};

#endif // FLUIDCELL_H
