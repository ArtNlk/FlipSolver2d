#include "benchruntable.h"

#include <OpenXLSX.hpp>

BenchRunTable::BenchRunTable()
{

}

void BenchRunTable::addStepTiming(const SolverStats &stats)
{
    m_currentSceneStats.push_back(stats);
}

void BenchRunTable::save()
{
    OpenXLSX::XLDocument doc;

    doc.create(m_outputFilePath.string(), true);

    for(auto& iter : m_allSceneStats)
    {
        const std::string& sceneName = iter.first;
        const StatVector& sceneStats = iter.second;

        doc.workbook().addWorksheet(sceneName);
        auto wks = doc.workbook().worksheet(sceneName);

        writeHeader(wks);

        int currentRow = 2;

        for(const SolverStats& stats : sceneStats)
        {
            writeRow(stats, wks, currentRow);
            currentRow++;
        }
    }

    doc.save();
    doc.close();
}

void BenchRunTable::finishScene(const std::string &sceneName)
{
    m_allSceneStats.insert({sceneName, m_currentSceneStats});
    m_currentSceneStats.clear();
}

void BenchRunTable::writeHeader(OpenXLSX::XLWorksheet &sheet)
{
    TableColumn column = static_cast<TableColumn>(0);
    do
    {
        sheet.cell(1,1+column).value() = getColumnHeader(column);
    }while(nextEnum<TableColumn,TABLE_COLUMN_COUNT>(column));
}

void BenchRunTable::writeRow(const SolverStats &stats, OpenXLSX::XLWorksheet &worksheet, int rowIdx)
{
    worksheet.cell(rowIdx,1+STEP_NUMBER_COLUMN).value() = rowIdx-1;
    worksheet.cell(rowIdx,1+SUBSTEP_COUNT_COLUMN).value() = stats.substepCount();
    worksheet.cell(rowIdx,1+TOTAL_FRAME_TIME_COLUMN).value() = stats.frameTime();
    worksheet.cell(rowIdx,1+ADVECTION_TIME_COLUMN).value() = stats.timings()[SolverStage::ADVECTION];
    worksheet.cell(rowIdx,1+DECOMPOSITION_TIME_COLUMN).value() = stats.timings()[SolverStage::DECOMPOSITION];
    worksheet.cell(rowIdx,1+DENSITY_TIME_COLUMN).value() = stats.timings()[SolverStage::DENSITY];
    worksheet.cell(rowIdx,1+PARTICLE_REBIN_TIME_COLUMN).value() = stats.timings()[SolverStage::PARTICLE_REBIN];
    worksheet.cell(rowIdx,1+PARTICLE_TO_GRID_TIME_COLUMN).value() = stats.timings()[SolverStage::PARTICLE_TO_GRID];
    worksheet.cell(rowIdx,1+GRID_UPDATE_TIME_COLUMN).value() = stats.timings()[SolverStage::GRID_UPDATE];
    worksheet.cell(rowIdx,1+AFTER_TRANSFER_TIME_COLUMN).value() = stats.timings()[SolverStage::AFTER_TRANSFER];
    worksheet.cell(rowIdx,1+PRESSURE_TIME_COLUMN).value() = stats.timings()[SolverStage::PRESSURE];
    worksheet.cell(rowIdx,1+VISCOSITY_TIME_COLUMN).value() = stats.timings()[SolverStage::VISCOSITY];
    worksheet.cell(rowIdx,1+REPRESSURE_TIME_COLUMN).value() = stats.timings()[SolverStage::REPRESSURE];
    worksheet.cell(rowIdx,1+PARTICLE_UPDATE_TIME_COLUMN).value() = stats.timings()[SolverStage::PARTICLE_UPDATE];
    worksheet.cell(rowIdx,1+PARTICLE_RESEED_TIME_COLUMN).value() = stats.timings()[SolverStage::PARTICLE_RESEED];
    worksheet.cell(rowIdx,1+PRESSURE_ITERS_COLUMN).value() = stats.pressureIterations();
    worksheet.cell(rowIdx,1+DENSITY_ITERS_COLUMN).value() = stats.densityIterations();
    worksheet.cell(rowIdx,1+VISCOSITY_ITERS_COLUMN).value() = stats.viscosityIterations();
}

std::string BenchRunTable::getColumnHeader(TableColumn column)
{
    switch(column)
    {
        case STEP_NUMBER_COLUMN:
            return "Step number";

        case SUBSTEP_COUNT_COLUMN:
            return "Substeps";

        case TOTAL_FRAME_TIME_COLUMN:
            return "Frame time";

        case ADVECTION_TIME_COLUMN:
            return "Advection";

        case DECOMPOSITION_TIME_COLUMN:
            return "Decomposition";

        case DENSITY_TIME_COLUMN:
            return "Density correction";

        case PARTICLE_REBIN_TIME_COLUMN:
            return "Particle rebin";

        case PARTICLE_TO_GRID_TIME_COLUMN:
            return "Particle to grid";

        case GRID_UPDATE_TIME_COLUMN:
            return "Grid update";

        case AFTER_TRANSFER_TIME_COLUMN:
            return "After transfer";

        case PRESSURE_TIME_COLUMN:
            return "Pressure";

        case VISCOSITY_TIME_COLUMN:
            return "Viscosity";

        case REPRESSURE_TIME_COLUMN:
            return "After-visc pressure";

        case PARTICLE_UPDATE_TIME_COLUMN:
            return "Particle update";

        case PARTICLE_RESEED_TIME_COLUMN:
            return "Particle reseeding";

        case PRESSURE_ITERS_COLUMN:
            return "Pressure iterations";

        case DENSITY_ITERS_COLUMN:
            return "Density iterations";

        case VISCOSITY_ITERS_COLUMN:
            return "Viscosity iterations";

        case TABLE_COLUMN_COUNT:
            return "INVALID COLUMN";
    }

    return "INVALID ENUM";
}
