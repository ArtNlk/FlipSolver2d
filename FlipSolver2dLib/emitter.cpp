#include "emitter.h"

Emitter::Emitter(float viscosity, Geometry2d geo):
    m_viscosity(viscosity),
    m_geometry(geo)
{

}

void Emitter::setViscosity(float viscosity)
{
    m_viscosity = viscosity;
}

float Emitter::viscosity()
{
    return m_viscosity;
}

Geometry2d &Emitter::geometry()
{
    return m_geometry;
}

void Emitter::setGeometry(Geometry2d& geometry)
{
    m_geometry = geometry;
}
