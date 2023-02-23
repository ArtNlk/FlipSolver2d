#ifndef EMITTER_H
#define EMITTER_H

#include "geometry2d.h"

class  Emitter
{
public:
    Emitter(Geometry2d& geo);

    void setViscosity(float viscosity);

    float viscosity() const;

    Geometry2d &geometry();

    void setGeometry(Geometry2d &newGeometry);

    float temperature() const;
    void setTemperature(float newTemperature);

    float concentrartion() const;
    void setConcentrartion(float newConcentrartion);

    float divergence() const;
    void setDivergence(float newDivergence);

    float fuel() const;
    void setFuel(float newFuel);

    Vertex velocity() const;
    void setVelocity(Vertex newVelocity);

    bool velocityTransfer() const;
    void setVelocityTransfer(bool velocityTransfer);

protected:
    float m_viscosity;
    float m_temperature;
    float m_concentrartion;
    float m_divergence;
    float m_fuel;
    bool m_transferVelocity;
    Vertex m_velocity;
    Geometry2d m_geometry;
};

#endif // EMITTER_H
