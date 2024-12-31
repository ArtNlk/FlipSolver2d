#ifndef BENCHMARKRUNNERAPPLICATION_H
#define BENCHMARKRUNNERAPPLICATION_H

#include <filesystem>

#include "benchruntable.h"

class BenchmarkRunnerApplication
{
public:
    BenchmarkRunnerApplication(int argc, char* argv[]);

    int run();

private:
    void runScene(std::filesystem::path scenePath);

    std::filesystem::path m_sceneInputPath;
    std::filesystem::path m_outputPath;

    BenchRunTable m_timingTable;

    int m_stepCount;
};

#endif // BENCHMARKRUNNERAPPLICATION_H
