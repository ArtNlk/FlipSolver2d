#ifndef SIMSETTINGS_H
#define SIMSETTINGS_H


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

protected:
    static SimSettings m_instance;

    double m_dt;
    double m_dx;
    double m_density;
};

#endif // SIMSETTINGS_H
