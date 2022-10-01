#include "emitter.h"

Emitter::Emitter(float viscosity, float temperature, float concentration, Geometry2d& geo):
    m_viscosity(viscosity),
    m_temperature(temperature),
    m_concentrartion(concentration),
    m_geometry(geo)
{

}

void Emitter::setViscosity(float viscosity)
{
    m_viscosity = viscosity;
}

float Emitter::viscosity() const
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

float Emitter::temperature() const
{
    return m_temperature;
}

void Emitter::setTemperature(float newTemperature)
{
    m_temperature = newTemperature;
}

float Emitter::concentrartion() const
{
    return m_concentrartion;
}

void Emitter::setConcentrartion(float newConcentrartion)
{
    m_concentrartion = newConcentrartion;
}
