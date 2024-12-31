#ifndef BENCHRUNTABLE_H
#define BENCHRUNTABLE_H

#include "XLSheet.hpp"
#include "flipsolver2d.h"

#include <filesystem>
#include <map>

class BenchRunTable
{
public:
    using StatVector = std::vector<SolverStats>;

    BenchRunTable();

    void setOutputFilePath(std::filesystem::path& outputPath)
    {
        m_outputFilePath = outputPath;
    }

    void addStepTiming(const SolverStats& stats);

    void save();

    void finishScene(const std::string& sceneName);

    enum TableColumn
    {
        STEP_NUMBER_COLUMN = 0,
        SUBSTEP_COUNT_COLUMN,
        TOTAL_FRAME_TIME_COLUMN,
        ADVECTION_TIME_COLUMN,
        DECOMPOSITION_TIME_COLUMN,
        DENSITY_TIME_COLUMN,
        PARTICLE_REBIN_TIME_COLUMN,
        PARTICLE_TO_GRID_TIME_COLUMN,
        GRID_UPDATE_TIME_COLUMN,
        AFTER_TRANSFER_TIME_COLUMN,
        PRESSURE_TIME_COLUMN,
        VISCOSITY_TIME_COLUMN,
        REPRESSURE_TIME_COLUMN,
        PARTICLE_UPDATE_TIME_COLUMN,
        PARTICLE_RESEED_TIME_COLUMN,
        PRESSURE_ITERS_COLUMN,
        DENSITY_ITERS_COLUMN,
        VISCOSITY_ITERS_COLUMN,
        TABLE_COLUMN_COUNT
    };

protected:
    void writeHeader(OpenXLSX::XLWorksheet& sheet);

    void writeRow(const SolverStats& stats, OpenXLSX::XLWorksheet& worksheet, int rowIdx);

    std::string getColumnHeader(TableColumn column);

    std::map<std::string,StatVector> m_allSceneStats;

    StatVector m_currentSceneStats;

    std::filesystem::path m_outputFilePath;
};

#endif // BENCHRUNTABLE_H
