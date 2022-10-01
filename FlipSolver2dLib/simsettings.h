#ifndef SIMSETTINGS_H
#define SIMSETTINGS_H

#include "geometry2d.h"

enum SimulationType : char {SIMULATION_FLUID, SIMULATION_GAS};

class SimSettings
{
public:
    SimSettings();

    static inline SimSettings& i()
    {
        return m_instance;
    }

    static inline float& stepDt()
    {
        return m_instance.m_stepDt;
    }

    static inline float& frameDt()
    {
        return m_instance.m_frameDt;
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

    static inline int& domainSizeI()
    {
        return m_instance.m_domainSizeI;
    }

    static inline int& domainSizeJ()
    {
        return m_instance.m_domainSizeJ;
    }

    static inline int& gridSizeI()
    {
        return m_instance.m_gridSizeI;
    }

    static inline int& gridSizeJ()
    {
        return m_instance.m_gridSizeJ;
    }

    static inline float& resolution()
    {
        return m_instance.m_resolution;
    }

    static inline int& fps()
    {
        return m_instance.m_fps;
    }

    static inline int& maxSubsteps()
    {
        return m_instance.m_maxSubsteps;
    }

    static inline float& picRatio()
    {
        return m_instance.m_picRatio;
    }

    static inline float& cflNumber()
    {
        return m_instance.m_cflNumber;
    }

    static inline SimulationType& simType()
    {
        return m_instance.m_simType;
    }

    static inline float& ambientTemp()
    {
        return m_instance.m_ambientTemperature;
    }

protected:
    static SimSettings m_instance;

    float m_stepDt;
    float m_frameDt;
    double m_dx;
    double m_density;
    unsigned int m_seed;
    int m_particlesPerCell;
    Vertex m_globalAcceleration;
    int m_domainSizeI;
    int m_domainSizeJ;
    int m_gridSizeI;
    int m_gridSizeJ;
    float m_resolution;
    int m_fps;
    int m_maxSubsteps;
    float m_picRatio;
    float m_cflNumber;
    float m_ambientTemperature;
    SimulationType m_simType;
};

#endif // SIMSETTINGS_H
