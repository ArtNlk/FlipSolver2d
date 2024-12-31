#include "benchmarkrunnerapplication.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "argumentparser.h"
#include "jsonscenereader.h"

namespace
{
static const char* ScenePathArgument = "-i";
static const char* OutputPathArgument = "-o";
static const char* StepCountArgument = "-s";
}

BenchmarkRunnerApplication::BenchmarkRunnerApplication(int argc, char *argv[])
{
    ArgumentParser args(argc, argv);

    m_sceneInputPath = args.getCmdOption(ScenePathArgument);
    m_outputPath = args.getCmdOption(OutputPathArgument);

    if(m_sceneInputPath.empty())
    {
        m_sceneInputPath = std::filesystem::current_path();
    }

    if(m_outputPath.empty())
    {
        m_outputPath = std::filesystem::current_path();
    }

    std::string stepCountStr = args.getCmdOption(StepCountArgument);

    if(stepCountStr.empty())
    {
        m_stepCount = 60;
    }
    else
    {
        m_stepCount = std::stoi(stepCountStr);
    }

    m_timingTable.setOutputFilePath(m_outputPath.append("Stats.xlsx"));
}

int BenchmarkRunnerApplication::run()
{
    if(!std::filesystem::is_directory(m_sceneInputPath) || m_sceneInputPath.extension() == ".json")
    {
        runScene(m_sceneInputPath);
        m_timingTable.save();
        return 0;
    }

    for(auto const& sceneFolderEntry : std::filesystem::directory_iterator(m_sceneInputPath))
    {
        if(std::filesystem::is_directory(sceneFolderEntry) || sceneFolderEntry.path().extension() != ".json")
        {
            continue;
        }

        runScene(sceneFolderEntry);
    }

    m_timingTable.save();

    return 0;
}

void BenchmarkRunnerApplication::runScene(std::filesystem::path scenePath)
{
    const std::string pathString = scenePath.string();

    std::shared_ptr solver = JsonSceneReader::loadJson(pathString);

    if(!solver)
    {
        std::cout << "Failed to load scene: " << pathString << '\n';
        return;
    }

    std::cout << "Starting scene: " << pathString << '\n';

    for(int step = 0; step < m_stepCount; step++)
    {
        solver->stepFrame();
        m_timingTable.addStepTiming(solver->timeStats());
    }

    std::cout << "Finished scene: " << pathString << '\n';

    m_timingTable.finishScene(scenePath.stem());
}
