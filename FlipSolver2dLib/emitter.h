#ifndef EMITTER_H
#define EMITTER_H

#include "geometry2d.h"

class Emitter
{
public:
    Emitter(float viscosity, Geometry2d geo);

    void setViscosity(float viscosity);

    float viscosity();

    Geometry2d &geometry();

    void setGeometry(Geometry2d &newGeometry);

protected:
    float m_viscosity;
    Geometry2d m_geometry;
};

#endif // EMITTER_H
