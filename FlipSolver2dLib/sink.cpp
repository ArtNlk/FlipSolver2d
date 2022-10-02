#include "sink.h"

Sink::Sink(float divergence, Geometry2d &geo) :
    m_geo(geo),
    m_divergence(divergence)
{

}

float Sink::divergence() const
{
    return m_divergence;
}

void Sink::setDivergence(float newDivergence)
{
    m_divergence = newDivergence;
}

Geometry2d &Sink::geo()
{
    return m_geo;
}

void Sink::setGeo(const Geometry2d &newGeo)
{
    m_geo = newGeo;
}
