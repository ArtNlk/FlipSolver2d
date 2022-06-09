#include "simsettings.h"

SimSettings SimSettings::m_instance;

SimSettings::SimSettings()
{
    m_domainSizeI = 50;
    m_domainSizeJ = 50;
    m_resolution = 50;
    if(m_domainSizeI > m_domainSizeJ)
    {
        m_dx = static_cast<float>(m_domainSizeI) / m_resolution;
        m_gridSizeI = m_resolution;
        m_gridSizeJ = (static_cast<float>(m_domainSizeJ) / static_cast<float>(m_domainSizeI)) * m_resolution;
    }
    else
    {
        m_dx = static_cast<float>(m_domainSizeJ) / m_resolution;
        m_gridSizeJ = m_resolution;
        m_gridSizeI = (static_cast<float>(m_domainSizeI) / static_cast<float>(m_domainSizeJ)) * m_resolution;
    }
    m_fps = 30;
    m_substeps = 1;
    m_density = 100;
    m_stepDt = 1.f/ (m_fps * m_substeps);
    m_frameDt = 1.f / m_fps;
    m_seed = 0;
    m_particlesPerCell = 4;
    m_globalAcceleration = Vertex(9.8f,0);
}
