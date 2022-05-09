#include "logger.h"

#include "dynamicuppertriangularsparsematrix.h"
#include "uppertriangularmatrix.h"

Logger::Logger()
{
    m_filePath = "./FluidSimLog.txt";
    m_logFileStream.open(m_filePath);
    if(m_logFileStream.fail())
    {
        std::cout << "Error opening log file!";
    }
}

Logger::~Logger()
{
    m_logFileStream.close();
}

Logger &operator<<(Logger &l, const DynamicUpperTriangularSparseMatrix &m)
{
    std::ofstream &s = l.stream();
    s << '[';
    for(int i = 0; i < m.sizeI(); i++)
    {
        for(int j = 0; j < m.sizeJ(); j++)
        {
            s << m.getValue(i,j) << ' ';
        }
        s << ";";
    }
    s << ']';

    return l;
}

Logger &operator<<(Logger &l, const UpperTriangularMatrix &m)
{
    std::ofstream &s = l.stream();
    s << '[';
    for(int i = 0; i < m.sizeI(); i++)
    {
        for(int j = 0; j < m.sizeJ(); j++)
        {
            s << m.getValue(i,j) << ' ';
        }
        s << ";";
    }
    s << ']';

    return l;
}

Logger &operator<<(Logger &l, const std::vector<double> &v)
{
    std::ofstream &s = l.stream();
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i] << ' ';
    }
    s << ']';

    return l;
}

Logger &operator<<(Logger &l, const std::vector<float> &v)
{
    std::ofstream &s = l.stream();
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i] << ' ';
    }
    s << ']';

    return l;
}

Logger &operator<<(Logger &l, const std::vector<int> &v)
{
    std::ofstream &s = l.stream();
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i] << ' ';
    }
    s << ']';

    return l;
}

Logger &operator<<(Logger &l, const std::string &str)
{
    l.stream() << str;
    return l;
}

Logger &operator<<(Logger &l, const double &v)
{
    l.stream() << v;
    return l;
}

Logger &operator<<(Logger &l, const float &v)
{
    l.stream() << v;
    return l;
}

Logger &operator<<(Logger &l, const int &v)
{
    l.stream() << v;
    return l;
}

Logger &operator<<(Logger &l, const char* str)
{
    l.stream() << str;
    return l;
}
