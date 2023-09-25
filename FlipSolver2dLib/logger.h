#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class DynamicMatrix;

class StaticMatrix;

class Logger
{
public:
    Logger();
    ~Logger();

    static Logger& instance()
    {
        static Logger instance;
        return instance;
    }

    std::ofstream &stream()
    {
        return m_logFileStream;
    }

    void close();

    Logger(Logger const&) = delete;
    void operator=(Logger const&) = delete;
protected:
    std::string m_filePath;
    std::ofstream m_logFileStream;
};

void binDump(StaticMatrix & m, std::string path);

Logger &operator<<(Logger &l, DynamicMatrix &m);

Logger &operator<<(Logger &l, const StaticMatrix &m);

Logger &operator<<(Logger &l, const std::vector<double> &v);

Logger &operator<<(Logger &l, const std::vector<float> &v);

Logger &operator<<(Logger &l, const std::vector<int> &v);

Logger &operator<<(Logger &l, const std::string &str);

Logger &operator<<(Logger &l, const double &v);

Logger &operator<<(Logger &l, const float &v);

Logger &operator<<(Logger &l, const int &v);

Logger &operator<<(Logger &l, const unsigned int &v);

Logger &operator<<(Logger &l, const char* str);

#define debug() Logger::instance().stream() << "\n"; Logger::instance()

#endif // LOGGER_H
