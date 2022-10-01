#ifndef EMITTER_H
#define EMITTER_H

#include "geometry2d.h"

class Emitter
{
public:
    Emitter(float viscosity, float temperature, float concentration, Geometry2d& geo);

    void setViscosity(float viscosity);

    float viscosity() const;

    Geometry2d &geometry();

    void setGeometry(Geometry2d &newGeometry);

    float temperature() const;
    void setTemperature(float newTemperature);

    float concentrartion() const;
    void setConcentrartion(float newConcentrartion);

protected:
    float m_viscosity;
    float m_temperature;
    float m_concentrartion;
    Geometry2d m_geometry;
};

#endif // EMITTER_H
