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

void Logger::close()
{
    std::cout << "closed";
    m_logFileStream.flush();
    m_logFileStream.close();
}

Logger &operator<<(Logger &l,DynamicUpperTriangularSparseMatrix &m)
{
#ifdef NUMPY_LOGGING
    std::ofstream &s = l.stream();
    s << '[';
    for(auto &sparseRow : m.data())
    {
        s << '[';
        for(auto &rowUnit : sparseRow)
        {
            s << '[' << rowUnit.first << ',' << rowUnit.second << "],";
        }
        s << "],";
    }
    s << ']';
#else
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
#endif

    s.flush();
    return l;
}

Logger &operator<<(Logger &l, const UpperTriangularMatrix &m)
{
    std::ofstream &s = l.stream();
#ifdef NUMPY_LOGGING
    s << '[';
    for(int i = 0; i < m.sizeI(); i++)
    {
        s << '[';
        for(int j = 0; j < m.sizeJ(); j++)
        {
            s << m.getValue(i,j);
            if(j != m.sizeJ() - 1) s << ',';
        }
        s << ']';
        if(i != m.sizeI() - 1) s << ',';
    }
    s << ']';
#else
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
#endif

    s.flush();
    return l;
}

Logger &operator<<(Logger &l, const std::vector<double> &v)
{
    std::ofstream &s = l.stream();
#ifdef NUMPY_LOGGING
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i];
        if(i != v.size() - 1) s << ',';
    }
    s << ']';
#else
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i] << ' ';
    }
    s << ']';
#endif
    s.flush();
    return l;
}

Logger &operator<<(Logger &l, const std::vector<float> &v)
{
    std::ofstream &s = l.stream();
#ifdef NUMPY_LOGGING
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i];
        if(i != v.size() - 1) s << ',';
    }
    s << ']';
#else
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i] << ' ';
    }
    s << ']';
#endif

    s.flush();
    return l;
}

Logger &operator<<(Logger &l, const std::vector<int> &v)
{
    std::ofstream &s = l.stream();
#ifdef NUMPY_LOGGING
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i];
        if(i != v.size() - 1) s << ',';
    }
    s << ']';
#else
    s << '[';
    for(int i = 0; i < v.size(); i++)
    {
        s << v[i] << ' ';
    }
    s << ']';
#endif

    s.flush();
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

Logger &operator<<(Logger &l, const unsigned int &v)
{
    l.stream() << v;
    return l;
}

Logger &operator<<(Logger &l, const char* str)
{
    l.stream() << str;
    return l;
}

void binDump(UpperTriangularMatrix &m, std::string path)
{
    std::ofstream s;
    s.open(path, std::ios_base::binary);
    for(int i = 0; i < m.sizeI(); i++)
    {
        for(int j = 0; j < m.sizeJ(); j++)
        {
            s << m.getValue(i,j);
        }
    }
    s.close();
}
