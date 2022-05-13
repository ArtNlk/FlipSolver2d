#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class DynamicUpperTriangularSparseMatrix;

class UpperTriangularMatrix;

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

    Logger(Logger const&) = delete;
    void operator=(Logger const&) = delete;
protected:
    std::string m_filePath;
    std::ofstream m_logFileStream;
};

Logger &operator<<(Logger &l, const DynamicUpperTriangularSparseMatrix &m);

Logger &operator<<(Logger &l, const UpperTriangularMatrix &m);

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
