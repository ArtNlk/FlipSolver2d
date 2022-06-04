#ifndef SIMSETTINGS_H
#define SIMSETTINGS_H

#include "geometry2d.h"

class SimSettings
{
public:
    SimSettings();

    static inline SimSettings& i()
    {
        return m_instance;
    }

    static inline double& dt()
    {
        return m_instance.m_dt;
    }

    static inline double& dx()
    {
        return m_instance.m_dx;
    }

    static inline double& density()
    {
        return m_instance.m_density;
    }

    static inline unsigned int& randomSeed()
    {
        return m_instance.m_seed;
    }

    static inline int& particlesPerCell()
    {
        return m_instance.m_particlesPerCell;
    }

    static inline Vertex& globalAcceleration()
    {
        return m_instance.m_globalAcceleration;
    }

protected:
    static SimSettings m_instance;

    double m_dt;
    double m_dx;
    double m_density;
    unsigned int m_seed;
    int m_particlesPerCell;
    Vertex m_globalAcceleration;
};

#endif // SIMSETTINGS_H
