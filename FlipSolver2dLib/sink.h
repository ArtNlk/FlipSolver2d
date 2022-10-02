#ifndef SINK_H
#define SINK_H

#include "geometry2d.h"

class Sink
{
public:
    Sink(float divergence, Geometry2d &geo);

    float divergence() const;
    void setDivergence(float newDivergence);

    Geometry2d &geo();
    void setGeo(const Geometry2d &newGeo);

protected:
    Geometry2d m_geo;

    float m_divergence;
};

#endif // SINK_H
