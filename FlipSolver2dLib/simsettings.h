#ifndef SIMSETTINGS_H
#define SIMSETTINGS_H


class SimSettings
{
public:
    SimSettings();

    inline SimSettings& i()
    {
        return m_instance;
    }

    inline double& dt()
    {
        return m_dt;
    }

    inline double& dx()
    {
        return m_dx;
    }

    inline double& density()
    {
        return m_density;
    }

protected:
    static SimSettings m_instance;

    double m_dt;
    double m_dx;
    double m_density;
};

#endif // SIMSETTINGS_H
