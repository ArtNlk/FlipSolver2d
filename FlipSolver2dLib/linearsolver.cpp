#include "linearsolver.h"

#include <algorithm>

#include <array>
#include <cmath>
#include <future>
#include <mutex>
#include <ostream>
#include <sstream>
#include <utility>
#include <vector>
#include <format>

#include "PressureIPPCoeficients.h"
#include "dynamicmatrix.h"
#include "grid2d.h"
#include "materialgrid.h"
#include "pressuredata.h"
#include "vmath.h"
#include "linearsolver_sse42.h"

#include <fstream>

#include "logger.h"

LinearSolver::LinearSolver(MaterialGrid &materialGrid, int maxMultigridDepth) :
    m_restrictionWeights({0}),
    m_prolongationWeights({0}),
    m_mainMaterialGrid(materialGrid)
{

    int sizeI = materialGrid.sizeI();
    int sizeJ = materialGrid.sizeJ();
    int subgridLevelCount = (std::min(sizeI,sizeJ) / 8) - 1;
    subgridLevelCount = std::min(maxMultigridDepth,subgridLevelCount);
    for(int level = 0; level < subgridLevelCount; level++)
    {
        sizeI /= 2;
        sizeJ /= 2;
        if(sizeI < 8 || sizeJ < 8)
        {
            break;
        }
        m_materialSubgrids.emplace_back(sizeI,sizeJ,FluidMaterial::SOLID);
        m_pressureGrids.emplace_back(sizeI,sizeJ,0.f,OOBStrategy::OOB_EXTEND);
        m_rhsGrids.emplace_back(sizeI,sizeJ,0.f,OOBStrategy::OOB_EXTEND);
        std::cout << "Added subgrid IxJ: " << sizeI << ' ' << sizeJ << std::endl;
    }
    const std::array<float, 4> weights = {1.f,3.f,3.f,1.f};
    int idx = 0;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            m_restrictionWeights[idx] = (weights[i] * weights[j]) / 64.f;
            m_prolongationWeights[idx] = (weights[i] * weights[j]) / 32.f;
            idx++;
        }
    }

}

int LinearSolver::solve(const IndexedPressureParameters &matrixIn,
                        const IndexedIPPCoefficients &precond,
                        std::vector<double> &result,
                        const std::vector<double> &vec,
                        int iterLimit,
                        double tol)
{
    result.assign(result.size(), 0);
    if (VOps::i().isZero(vec)) {
        std::cout << "Solver skipped zeros vector" << '\n';
        return true;
    }
    updateSubgrids();

    //debug() << "mat=" << matrix;
    std::vector<double> residual(vec);
    std::vector<double> aux = residual;
    //applyICPrecond(precond,residual,aux);
    //applyIPPrecond(&matrix,&residual,&aux);
    //precond->apply(residual,aux,data);
    std::vector<double> search = aux;
    double sigma = VOps::i().dot(aux, residual);
    double err = 0.0;
    for (int i = 0; i < iterLimit; i++)
    {
        matrixIn.multiply(search, aux);
        double alpha = sigma / (VOps::i().dot(aux, search) + 1e-8);
        VOps::i().addMul(result, result, search, alpha);
        VOps::i().subMul(residual, residual, aux, alpha);
        err = VOps::i().maxAbs(residual);

        if (err <= tol)
        {
            debug() << "[SOLVER] Solver done, iter = " << i << " err = " << err;
            std::cout << "Solver done, iter = " << i << " err = " << err << '\n';

            return i;
        }

        //aux = residual;
        applyMGPrecond(residual, aux);
        //precond.multiply(residual, aux);
        //applyICPrecond(precond,residual,aux);
        //applyIPPrecond(&matrix,&residual,&aux);
        double newSigma = VOps::i().dot(aux, residual);
        double beta = newSigma / (sigma);
        VOps::i().addMul(search, aux, search, beta);
        sigma = newSigma;
    }

    debug() << "Solver iter exhaustion, err = " << err;
    std::cout << "Solver iter exhaustion, err = " << err << '\n';
    return iterLimit;
}

void LinearSolver::applyMGPrecond(const std::vector<double> &in, std::vector<double> &out)
{
    vCycle(out,in);
}

void LinearSolver::updateSubgrids()
{
    // std::ofstream file("grids.txt", std::ios::out);

    // file << m_mainMaterialGrid.toString() << '\n';
    // for(int i = 0; i < m_materialSubgrids.size(); i++)
    // {
    //     file << m_materialSubgrids[i].toString() << '\n';
    // }

    propagateMaterialGrid(m_mainMaterialGrid, m_materialSubgrids[0]);
    for(int i = 0; i < (m_materialSubgrids.size() - 1); i++)
    {
        propagateMaterialGrid(m_materialSubgrids[i], m_materialSubgrids[i+1]);
    }

    // file << m_mainMaterialGrid.toString() << '\n';
    // for(int i = 0; i < m_materialSubgrids.size(); i++)
    // {
    //     file << m_materialSubgrids[i].toString() << '\n';
    // }
    // file.close();
}

void LinearSolver::propagateMaterialGrid(const MaterialGrid &fineGrid, MaterialGrid &coarseGrid)
{
    //neumann - solid
    //dirichlet - air
    //interior - fluid
    const std::vector<FluidMaterial> &fineGridData = fineGrid.data();
    const int gridDataSize = fineGridData.size();
    for(int i = 0; i < coarseGrid.sizeI(); i++)
    {
        for(int j = 0; j < coarseGrid.sizeJ(); j++)
        {
            std::array<int,4> childCellIdxs = getFineGridChildIdxs(fineGrid,i,j);
            FluidMaterial m = FluidMaterial::SOLID;
            for(int idx : childCellIdxs)
            {
                if(idx < 0 || idx >= gridDataSize)
                {
                    continue;
                }
                if(fineGridData[idx] == FluidMaterial::EMPTY)
                {
                    m = FluidMaterial::EMPTY;
                    break;
                }
                if(fineGridData[idx] == FluidMaterial::FLUID)
                {
                    m = FluidMaterial::FLUID;
                }
            }
            coarseGrid.at(i,j) = m;
        }
    }
}
void LinearSolver::restrictGrid(const MaterialGrid& fineMaterials, const MaterialGrid &coarseMaterials, const Grid2d<double> &fineGrid, Grid2d<double> coarseGrid)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(coarseMaterials.data().size(),128);
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::restrictGridThread,this,range,fineMaterials,coarseMaterials,
                                 std::ref(fineGrid),std::ref(coarseGrid));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    ThreadPool::i()->wait();
}
void LinearSolver::restrictGridThread(const Range range, const MaterialGrid& fineMaterials, const MaterialGrid &coarseMaterials, const Grid2d<double> &fineGrid, Grid2d<double> coarseGrid)
{
    const std::vector<double> &fineGridData = fineGrid.data();
    const int gridDataSize = fineGridData.size();
    for(int idx = range.start; idx < range.end; idx++)
    {
        Index2d i2d = coarseMaterials.index2d(idx);
        std::array<int,16> stencilCellIdxs = getFineGridStencilIdxs(fineGrid,i2d.i,i2d.j);
        coarseGrid.at(i2d.i,i2d.j) = 0.f;
        if(!coarseMaterials.isFluid(i2d.i,i2d.j))
        {
            continue;
        }
        int weightIdx = 0;
        for(int idx : stencilCellIdxs)
        {
            if(idx < 0 || idx >= gridDataSize || fineMaterials.data()[idx] != FluidMaterial::FLUID)
            {
                weightIdx++;
                continue;
            }
            coarseGrid.at(i2d.i,i2d.j) += fineGridData[idx] * m_restrictionWeights[weightIdx];
            weightIdx++;
        }
    }
}
void LinearSolver::prolongateGrid(const MaterialGrid &coarseMaterials, const MaterialGrid &fineMaterials, const Grid2d<double> &coarseGrid, std::vector<double>& fineGridData)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(fineMaterials.data().size(),128);
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::prolongateGridThread,this,range,coarseMaterials,fineMaterials,
                                 std::ref(coarseGrid),std::ref(fineGridData));
    }
    ThreadPool::i()->wait();
}
void LinearSolver::prolongateGridThread(const Range range, const MaterialGrid &coarseMaterials, const MaterialGrid &fineMaterials, const Grid2d<double> &coarseGrid, std::vector<double> &fineGridData)
{
    const std::vector<double> &coarseGridData = coarseGrid.data();
    const int gridDataSize = coarseGridData.size();
    for(int idx = range.start; idx < range.end; idx++)
    {
        Index2d i2d = fineMaterials.index2d(idx);
        if(!fineMaterials.isFluid(i2d))
        {
            continue;
        }
        std::array<int, 4> coarseNeghbors = getCoarseProlongIdxs(coarseGrid,i2d.i,i2d.j);
        std::array<float,4> coarseWeights = getProlongationWeights(i2d.i,i2d.j);
        for(int idx = 0; idx < 4; idx++)
        {
            if(idx < 0 || idx >= gridDataSize || coarseMaterials.data()[idx] != FluidMaterial::FLUID)
            {
                continue;
            }
            fineGridData[fineMaterials.linearIndex(i2d.i,i2d.j)] +=
                coarseGridData[coarseNeghbors[idx]] * coarseWeights[idx];
        }
    }
}
void LinearSolver::premaskPressures(const MaterialGrid &materials, std::vector<double> &pressures)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(pressures.size(),128);
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(LinearSolver::premaskPressuresThread,range,
                                 std::ref(materials),std::ref(pressures));
    }
    ThreadPool::i()->wait();
}
void LinearSolver::premaskPressuresThread(const Range range, const MaterialGrid &materials, std::vector<double> &pressures)
{
    for(int i = range.start; i < range.end; i++)
    {
        //pressures[i] *= !solidTest(materials.data()[i]);
        pressures[i] = !fluidTest(materials.data()[i]) ? 0.0 : pressures[i];
    }
}
void LinearSolver::dampedJacobi(const MaterialGrid &materials, std::vector<double> &pressures, const std::vector<double> &rhs)
{
    std::vector<double> temp = pressures;
    std::vector<Range> ranges = ThreadPool::i()->splitRange(pressures.size(),128);
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::dampedJacobiThread,this,range,std::ref(materials),
                                 std::ref(temp),std::ref(pressures),std::ref(rhs));
    }
    ThreadPool::i()->wait();
    std::swap(temp,pressures);
    //    for(int i = 0; i < pressures.size(); i++)
    //    {
    //        pressures[i] += temp[i] * tune;
    //    }
    //    for(Range& range : ranges)
    //    {
    //        ThreadPool::i()->enqueue(&LinearSolver::vaddmul,this,range,std::ref(pressures),
    //                                 std::ref(temp),tune);
    //    }
    //    ThreadPool::i()->wait();
}

void LinearSolver::dampedJacobiThread(LinearSolver* solver, const Range range, const MaterialGrid &materials, std::vector<double> &vout, const std::vector<double> &pressures, const std::vector<double> &rhs)
{
    const double tune = 0.66666;
    const size_t dataSize = materials.data().size();
    for(int i = range.start; i < range.end; i++)
    {
        if(!materials.isFluid(static_cast<ssize_t>(i)))
        {
            vout[i] += (rhs[i]) * tune;
            continue;
        }
        //const auto weights = solver->getMultigridMatrixEntriesForCell(materials,i);
        const auto neighbors = materials.immidiateNeighbors(i);
        double result = 0.0;
        double diag = 0.0;
        for(int idx : neighbors)
        {
            if(idx < 0 || idx >= dataSize)
            {
                continue;
            }
            result -= pressures[idx]*materials.isFluid(static_cast<ssize_t>(idx));
            diag += !materials.isSolid(static_cast<ssize_t>(idx));
        }

        vout[i] += ((rhs[i]-result)/diag) * tune;
    }
}

void LinearSolver::vCycle(std::vector<double> &vout, const std::vector<double> &vin)
{
    const int strokeSmoothIter = 1;
    vout.assign(vin.size(),0.f);
    const int subLevelCount = m_materialSubgrids.size();
    //Down stroke
    for(int level = -1; level < static_cast<int>(m_materialSubgrids.size() - 1); level++)
    {
        if(level != -1)
        {
            Grid2d<double> temp(m_materialSubgrids[level].sizeI(),m_materialSubgrids[level].sizeJ(),0.f);
            //premaskPressures(m_materialSubgrids[level], m_pressureGrids[level].data());
            for(int i = 0; i < strokeSmoothIter; i++)
            {
                dampedJacobi(m_materialSubgrids[level],m_pressureGrids[level].data(),m_rhsGrids[level].data());
            }
            borderSmooth(m_materialSubgrids[level],m_materialSubgrids[level+1],m_rhsGrids[level].data(),m_pressureGrids[level].data(), 2<<level);
            multigridSubMatmul(m_materialSubgrids[level],
                               m_rhsGrids[level].data(),
                               m_pressureGrids[level].data(),
                               temp.data());
            restrictGrid(m_materialSubgrids[level],m_materialSubgrids[level+1],temp,m_rhsGrids[level+1]);
            m_pressureGrids[level+1].data().assign(m_pressureGrids[level+1].data().size(),0.f);
        }
        else
        {
            Grid2d<double> temp(m_mainMaterialGrid.sizeI(),m_mainMaterialGrid.sizeJ(),0.f);
            //premaskPressures(m_mainMaterialGrid, vout);

            for(int i = 0; i < strokeSmoothIter; i++)
            {
                dampedJacobi(m_mainMaterialGrid,vout,vin);
            }
            borderSmooth(m_mainMaterialGrid,m_materialSubgrids[0],vin,vout,2);
            multigridSubMatmul(m_mainMaterialGrid,
                               vin,
                               vout,
                               temp.data());
            restrictGrid(m_mainMaterialGrid,m_materialSubgrids[0],temp,m_rhsGrids[0]);
            m_pressureGrids[0].data().assign(m_pressureGrids[0].data().size(),0.f);
        }
    }

    //Solve coarse
    for(int i = 0; i < 50; i++)
    {
        //premaskPressures(m_materialSubgrids[subLevelCount-1], m_pressureGrids[subLevelCount-1].data());
        dampedJacobi(m_materialSubgrids[subLevelCount-1],
                     m_pressureGrids[subLevelCount-1].data(),
                     m_rhsGrids[subLevelCount-1].data());
    }
    //Up stroke
    for(int level = m_materialSubgrids.size() - 2; level >= -1; level--)
    {
        if(level != -1)
        {
            prolongateGrid(m_materialSubgrids[level+1],m_materialSubgrids[level],m_pressureGrids[level+1],m_pressureGrids[level].data());
            //premaskPressures(m_materialSubgrids[level], m_pressureGrids[level].data());
            borderSmoothReverse(m_materialSubgrids[level],m_materialSubgrids[level+1],m_rhsGrids[level].data(),m_pressureGrids[level].data(), 2<<level);
            for(int i = 0; i < strokeSmoothIter; i++)
            {
            dampedJacobi(m_materialSubgrids[level],m_pressureGrids[level].data(),m_rhsGrids[level].data());
            }
        }
        else
        {
            prolongateGrid(m_materialSubgrids[0],m_mainMaterialGrid,m_pressureGrids[0],vout);
            //premaskPressures(m_mainMaterialGrid, vout);
            borderSmoothReverse(m_mainMaterialGrid,m_materialSubgrids[0],vin,vout,2);
            for(int i = 0; i < strokeSmoothIter; i++)
            {
                dampedJacobi(m_mainMaterialGrid,vout,vin);
            }
        }
    }
}

void LinearSolver::multigridSubMatmul(const MaterialGrid &materials, const std::vector<double> &vsub, const std::vector<double> &vmul, std::vector<double> &vout)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(vout.size(),8*8);
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::multigridSubMatmulThread,this,range,std::ref(materials),
                                 std::ref(vsub),std::ref(vmul),std::ref(vout));
    }
    ThreadPool::i()->wait();
}

void LinearSolver::multigridSubMatmulThread(const Range range, const MaterialGrid &materials, const std::vector<double> &vsub, const std::vector<double> &vmul, std::vector<double> &vout)
{
    const size_t dataSize = materials.data().size();
    for(int idx = range.start; idx < range.end; idx++)
    {
        //const auto weights = solver->getMultigridMatrixEntriesForCell(materials,i);
        const auto neighbors = materials.immidiateNeighbors(idx);
        double result = 0.0;
        if(materials.isFluid(static_cast<ssize_t>(idx)))
        {
            double diag = 0.0;
            for(int neighborIdx : neighbors)
            {
                if(neighborIdx < 0 || neighborIdx >= dataSize)
                {
                    continue;
                }
                result -= vmul[neighborIdx]*materials.isFluid(static_cast<ssize_t>(neighborIdx));
                diag += !materials.isSolid(static_cast<ssize_t>(neighborIdx));
            }
            result += diag * vmul[idx];
        }

        vout[idx] = vsub[idx] - result;
    }
}

void LinearSolver::borderSmooth(const MaterialGrid &currentMaterialGrid,
                                const LinearIndexable2d &coarserIndexer,
                                const std::vector<double> &rhsData,
                                std::vector<double> &smoothedPressureData, int iterCount)
{
    return;
    Grid2d<bool> maskGrid = generateBorderMask(currentMaterialGrid, coarserIndexer);

    std::ofstream file(std::format("maskGrid_{}x{}.txt",maskGrid.sizeI(), maskGrid.sizeJ()), std::ios::out);

    file << maskGrid.toString() << '\n';

    file.close();

    for(int smoothIter = 0; smoothIter < iterCount; smoothIter++)
    {
        for(ssize_t idx = 0; idx < maskGrid.linearSize(); idx++)
        {
            if(!maskGrid.data().at(idx) || !currentMaterialGrid.isFluid(static_cast<ssize_t>(idx)))
            {
                smoothedPressureData[idx] += (rhsData[idx]);
                continue;
            }
            //const auto weights = solver->getMultigridMatrixEntriesForCell(materials,i);
            const auto neighbors = currentMaterialGrid.immidiateNeighbors(idx);
            double result = 0.0;
            double diag = 0.0;
            for(int neighborIdx : neighbors)
            {
                if(neighborIdx < 0 || neighborIdx >= maskGrid.data().size())
                {
                    continue;
                }
                result -= smoothedPressureData[neighborIdx]*currentMaterialGrid.isFluid(static_cast<ssize_t>(neighborIdx));
                diag += !currentMaterialGrid.isSolid(static_cast<ssize_t>(neighborIdx));
            }

            smoothedPressureData[idx] += ((rhsData[idx]-result)/diag);
        }
    }
}

void LinearSolver::borderSmoothReverse(const MaterialGrid &currentMaterialGrid,
                                       const LinearIndexable2d &coarserIndexer,
                                       const std::vector<double> &rhsData,
                                       std::vector<double> &smoothedPressureData,
                                       int iterCount)
{
    return;
    Grid2d<bool> maskGrid = generateBorderMask(currentMaterialGrid, coarserIndexer);

    std::ofstream file(std::format("upstr_maskGrid_{}x{}.txt",maskGrid.sizeI(), maskGrid.sizeJ()), std::ios::out);

    file << maskGrid.toString() << '\n';

    file.close();

    for(int smoothIter = 0; smoothIter < iterCount; smoothIter++)
    {
        for(ssize_t idx = maskGrid.linearSize()-1; idx >= 0; idx--)
        {
            if(!maskGrid.data().at(idx) || !currentMaterialGrid.isFluid(static_cast<ssize_t>(idx)))
            {
                smoothedPressureData[idx] += (rhsData[idx]);
                continue;
            }
            //const auto weights = solver->getMultigridMatrixEntriesForCell(materials,i);
            const auto neighbors = currentMaterialGrid.immidiateNeighbors(idx);
            double result = 0.0;
            double diag = 0.0;
            for(int neighborIdx : neighbors)
            {
                if(neighborIdx < 0 || neighborIdx >= maskGrid.data().size())
                {
                    continue;
                }
                result -= smoothedPressureData[neighborIdx]*currentMaterialGrid.isFluid(static_cast<ssize_t>(neighborIdx));
                diag += !currentMaterialGrid.isSolid(static_cast<ssize_t>(neighborIdx));
            }

            smoothedPressureData[idx] += ((rhsData[idx]-result)/diag);
        }
    }
}

Grid2d<bool> LinearSolver::generateBorderMask(const MaterialGrid &currentMaterialGrid,
                                              const LinearIndexable2d &coarserIndexer)
{
    Grid2d<bool> output(currentMaterialGrid.sizeI(), currentMaterialGrid.sizeJ(), false);

    for(size_t maskIdx = 0; maskIdx < currentMaterialGrid.linearSize(); maskIdx++)
    {
        if(!currentMaterialGrid.isFluid(static_cast<ssize_t>(maskIdx)))
        {
            continue;
        }

        Index2d i2d = output.index2d(maskIdx);
        std::array<int, 16> checkedIdxs = getBorderCheckIndexes(coarserIndexer, currentMaterialGrid, i2d.i, i2d.j);

        for(int checkIdx : checkedIdxs)
        {
            if(checkIdx == -1)
            {
                break;
            }

            if(!currentMaterialGrid.isFluid(static_cast<ssize_t>(checkIdx)))
            {
                output.data()[maskIdx] = true;
                break;
            }
        }
    }

    return output;
}

std::array<int, 16> LinearSolver::getBorderCheckIndexes(const LinearIndexable2d &coarserGridIndexer,
                                                        const LinearIndexable2d &currentGridIndexer,
                                                        int i,
                                                        int j)
{
    std::array<int, 16> output{-1};

    std::array<int,4> coarseStencilIndexes = getCoarseProlongIdxs(coarserGridIndexer, i, j);

    size_t outputIndex = 0;

    for(int idx : coarseStencilIndexes)
    {
        if(idx < 0 || idx >= coarserGridIndexer.linearSize())
        {
            continue;
        }

        Index2d i2d = coarserGridIndexer.index2d(idx);
        std::array<int, 4> fineChildren = getFineGridChildIdxs(currentGridIndexer, i2d.i, i2d.j);

        for(int fineIdx : fineChildren)
        {
            if(fineIdx < 0 || fineIdx >= currentGridIndexer.linearSize())
            {
                continue;
            }

            output[outputIndex] = fineIdx;
            outputIndex++;
        }
    }

    return output;
}

std::pair<std::array<int,4>,std::array<double,4>> LinearSolver::getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int i, int j)
{
    return getMultigridMatrixEntriesForCell(materials, materials.linearIndex(i,j));
}

std::pair<std::array<int,4>,std::array<double,4>> LinearSolver::getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int linearIdx)
{
    std::pair<std::array<int,4>,std::array<double,4>> output;
    std::array<ssize_t,4> neighbors = materials.immidiateNeighbors(linearIdx);
    const std::vector<FluidMaterial> &materialsData = materials.data();
    const auto dataSize = materialsData.size();
    int outputIdx = 0;
    for(ssize_t idx : neighbors)
    {
        if(idx < 0 || idx >= dataSize)
        {
            output.first[outputIdx] = 0;
            output.second[outputIdx] = 0.;
            outputIdx++;
            continue;
        }
        const double weight = solidTest(materialsData[idx]) ? 0. : 1.;
        output.first[outputIdx] = idx;
        output.second[outputIdx] = weight;
        outputIdx++;
    }
    return output;
}
std::array<int, 16> LinearSolver::getFineGridStencilIdxs(const LinearIndexable2d &fineGridIndexer,int iCoarse, int jCoarse)
{
    std::array<int, 16> output{0};
    Index2d topLeftCorner(iCoarse*2 - 1, jCoarse*2 - 1);
    int topLeftLinearIdx = fineGridIndexer.linearIndex(topLeftCorner);
    int idx = 0;
    for(int iOffset = 0; iOffset < 4; iOffset++)
    {
        for(int jOffset = 0; jOffset < 4; jOffset++)
        {
            output[idx] = fineGridIndexer.linearIdxOfOffset(topLeftLinearIdx,iOffset, jOffset);
            idx++;
        }
    }
    return output;
}
std::array<int, 4> LinearSolver::getFineGridChildIdxs(const LinearIndexable2d &fineGridIndexer, int iCoarse, int jCoarse)
{
    std::array<int, 4> output{0};
    Index2d topLeftCorner(iCoarse*2, jCoarse*2);
    output[0] = fineGridIndexer.linearIndex(topLeftCorner);
    output[1] = fineGridIndexer.linearIdxOfOffset(output[0],0, 1);
    output[2] = fineGridIndexer.linearIdxOfOffset(output[0],1, 0);
    output[3] = fineGridIndexer.linearIdxOfOffset(output[0],1, 1);
    return output;
}

std::array<float, 4> LinearSolver::getProlongationWeights(int finei, int finej)
{
    int coarsePos = finei % 2;
    coarsePos |= (finej % 2) << 1;
    switch(coarsePos)
    {
    case 0:
        return {m_prolongationWeights[5],
                m_prolongationWeights[15],
                m_prolongationWeights[13],
                m_prolongationWeights[7]};
        break;
    case 1:
        return {m_prolongationWeights[6],
                m_prolongationWeights[12],
                m_prolongationWeights[14],
                m_prolongationWeights[4]};
        break;
    case 2:
        return {m_prolongationWeights[9],
                m_prolongationWeights[3],
                m_prolongationWeights[1],
                m_prolongationWeights[11]};
        break;
    case 3:
        return {m_prolongationWeights[10],
                m_prolongationWeights[0],
                m_prolongationWeights[2],
                m_prolongationWeights[8]};
        break;
    }
    return {0,0,0,0};
}
std::array<int, 4> LinearSolver::getCoarseProlongIdxs(const LinearIndexable2d &coarseGridIndexer, int finei, int finej)
{
    std::array<int,4> output;
    int iOffset = finei % 2 == 0 ? -1 : 1;
    int jOffset = finej % 2 == 0 ? -1 : 1;
    int currentCellIdx = coarseGridIndexer.linearIndex(finei/2,finej/2);
    output[0] = currentCellIdx;
    output[1] = coarseGridIndexer.linearIdxOfOffset(currentCellIdx,iOffset,jOffset);
    output[2] = coarseGridIndexer.linearIdxOfOffset(currentCellIdx,iOffset,0);
    output[3] = coarseGridIndexer.linearIdxOfOffset(currentCellIdx,0,jOffset);
    return output;
}
