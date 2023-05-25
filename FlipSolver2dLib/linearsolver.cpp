#include "linearsolver.h"

#include <algorithm>

#include <array>
#include <cmath>
#include <mutex>
#include <sstream>
#include <utility>
#include <vector>

#include "dynamicuppertriangularsparsematrix.h"
#include "grid2d.h"
#include "materialgrid.h"
#include "simddispatcher.h"
#include "vmath.h"
#include "linearsolver_sse42.h"

#include "logger.h"

const double LinearSolver::m_tol = 1.0e-6;

LinearSolver::LinearSolver(MaterialGrid &materialGrid, int maxMultigridDepth) :
    SimdDispatcher(),
    m_restrictionWeights({0}),
    m_prolongationWeights({0}),
    m_mainMaterialGrid(materialGrid),
    m_dampedJacobiThread(&LinearSolver::dampedJacobiThread)
{
    switch(m_simdLevel)
    {
    case SIMD_LEVEL_NONE:
        break;
    case SIMD_LEVEL_SSE42:
        m_dampedJacobiThread = &LinearSolver_sse42::dampedJacobiThread;
        break;
    case SIMD_LEVEL_SSE4a_XOP_FMA:
    case SIMD_LEVEL_AVX:
    case SIMD_LEVEL_AVX2:
    case SIMD_LEVEL_AVX512:
        break;
    }

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

bool LinearSolver::solve(const DynamicUpperTriangularSparseMatrix &matrixIn, std::vector<double> &result, const std::vector<double> &vec, int iterLimit)
{
    result.assign(result.size(),0);
    if (VOps::i().isZero(vec))
    {
        return true;
    }
    DynamicUpperTriangularSparseMatrix precond = calcPrecond(matrixIn);
    UpperTriangularMatrix matrix(matrixIn);

    //debug() << "mat=" << matrix;
    std::vector<double> residual = vec;
    std::vector<double> aux = vec;
    applyICPrecond(precond,residual,aux);
    std::vector<double> search = aux;
    double sigma = VOps::i().dot(aux,residual);
    double err = 0.0;
    for (int i = 0; i < iterLimit; i++)
    {
        aux = matrix * search;
        double alpha = sigma/(VOps::i().dot(aux,search));
        VOps::i().addMul(result,result,search,alpha);
        VOps::i().subMul(residual,residual,aux,alpha);
        err = VOps::i().maxAbs(residual);
//        if(i % 5 == 0)
//        {
//            std::cout << "Solver: " << i << " : " << err << "\n";
//            debug() << "Solver: " << i << " : " << err;
//        }
        if (err <= m_tol)
        {
            debug() << "[SOLVER] Solver done, iter = " << i << " err = " << err;
            std::cout << "Solver done, iter = " << i << " err = " << err << '\n';
            return true;
        }
        aux = residual;
        applyICPrecond(precond,residual,aux);
        double newSigma = VOps::i().dot(aux,residual);
        double beta = newSigma/(sigma);
        VOps::i().addMul(search,aux,search,beta);
        sigma = newSigma;
    }

    debug() << "Solver iter exhaustion, err = " << err;
    std::cout << "Solver iter exhaustion, err = " << err << '\n';
    return false;
}

bool LinearSolver::mfcgSolve(MatElementProvider elementProvider, std::vector<double> &result, const std::vector<double> &vec, int iterLimit)
{
    result.assign(result.size(),0);
    if (VOps::i().isZero(vec))
    {
        return true;
    }
    updateSubgrids();

    //vec - b
    std::vector<double> residual = vec; //r
    std::vector<double> aux = vec; //z
    applyMGPrecond(residual,aux);
    std::vector<double> search = aux; //p
    double sigma = VOps::i().dot(aux,residual);
    double err = 0.0;
    for (int i = 0; i < iterLimit; i++)
    {
        nomatVMul(elementProvider,search,aux);
        double alpha = sigma/(VOps::i().dot(aux,search));
        VOps::i().subMul(residual,residual,aux,alpha);
        err = VOps::i().maxAbs(residual);
        if (err <= m_tol)
        {
            //debug() << "[SOLVER] Solver done, iter = " << i << " err = " << err;
            std::cout << "MFSolver done, iter = " << i << " err = " << err << '\n';
            return true;
        }
        //aux = residual;
        applyMGPrecond(residual,aux);
        double newSigma = VOps::i().dot(aux,residual);
        //std::cout << "New sigma:" << newSigma << std::endl;
        double beta = newSigma/(sigma);
        VOps::i().addMul(result,result,search,alpha);
        VOps::i().addMul(search,aux,search,beta);
        sigma = newSigma;
    }

    //debug() << "Solver iter exhaustion, err = " << err;
    //std::cout << "Solver iter exhaustion, err = " << err << '\n';
    return false;
}

void LinearSolver::nomatVMul(MatElementProvider elementProvider, const std::vector<double> &vin, std::vector<double> &vout)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(vin.size());

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::nomatVMulThread,this,range,elementProvider,
                          std::ref(vin),std::ref(vout));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    ThreadPool::i()->wait();

//    for(int rowIdx = 0; rowIdx < vin.size(); rowIdx++)
//    {
//        SparseMatRowElements neighbors = elementProvider(rowIdx);
//        double mulSum = 0.0;
//        for(auto &matElement : neighbors)
//        {
//            if(matElement.first >= 0 && matElement.first < vin.size())
//            {
//                mulSum += vin.at(matElement.first) * matElement.second;
//            }
//        }
//        vout[rowIdx] = mulSum;
//    }
}

void LinearSolver::nomatVMulThread(Range range, MatElementProvider elementProvider, const std::vector<double> &vin, std::vector<double> &vout)
{
    for(unsigned int rowIdx = range.start; rowIdx < range.end; rowIdx++)
    {
        SparseMatRowElements neighbors = elementProvider(rowIdx);
        double mulSum = 0.0;
        for(auto &matElement : neighbors)
        {
            if(matElement.first >= 0 && matElement.first < vin.size())
            {
                mulSum += vin[matElement.first] * matElement.second;
            }
        }
        vout[rowIdx] = mulSum;
    }
}

void LinearSolver::applyMGPrecond(const std::vector<double> &in, std::vector<double> &out)
{
    vCycle(out,in);
}

void LinearSolver::applyICPrecond(const DynamicUpperTriangularSparseMatrix &precond, const std::vector<double> &in, std::vector<double> &out)
{
    out = in;

    std::vector<SparseRow> &rows = const_cast<DynamicUpperTriangularSparseMatrix&>(precond).data();

    for(int i = 0; i < out.size(); i++)
    {
        if(std::abs(precond.getValue(i,i)) > 1e-4)
        {
            auto& currentRow = rows[i];
            out[i] /= precond.getValue(i,i);

            for(int elementIdx = 1; elementIdx < currentRow.size(); elementIdx++)
            {
                int j = currentRow[elementIdx].first;
                if(j > i)
                {
                    out[j] = out[j] - precond.getValue(i,j) * out[i];
                }
            }
        }
    }

    for(int i = out.size() - 1; i >= 0; i--)
    {
        if(std::abs(precond.getValue(i,i)) > 1e-4)
        {
            auto& currentRow = rows[i];
            for(int elementIdx = 1; elementIdx < currentRow.size(); elementIdx++)
            {
                int j = currentRow[elementIdx].first;
                if(j>i)
                {
                    out[i] = out[i] - precond.getValue(i,j) * out[j];
                }
            }
            out[i] /= precond.getValue(i,i);
        }
    }
}

DynamicUpperTriangularSparseMatrix LinearSolver::calcPrecond(const DynamicUpperTriangularSparseMatrix &matrix)
{
    float tuning = 0.97f;
    float safety = 0.25f;
    int n = matrix.size();
    DynamicUpperTriangularSparseMatrix output(n);
    matrix.copyUpperTriangleTo(output);

    std::vector<SparseRow> &rowArray = output.data();

    for(int rowIdx = 0; rowIdx < rowArray.size(); rowIdx++)
    {
        if(rowArray[rowIdx][0].second != 0)
        {
            auto& currentRow = rowArray[rowIdx];
            if(output.getValue(rowIdx,rowIdx) < matrix.getValue(rowIdx,rowIdx) * safety)
            {
                output.setValue(rowIdx,rowIdx,sqrt(matrix.getValue(rowIdx,rowIdx)));
            }
            else
            {
                output.setValue(rowIdx,rowIdx,sqrt(output.getValue(rowIdx,rowIdx)));
            }

            for(int elementIdx = 1; elementIdx < currentRow.size(); elementIdx++)
            {
                currentRow[elementIdx].second = currentRow[elementIdx].second / output.getValue(rowIdx,rowIdx);
            }

            for(int elementIdx = 1; elementIdx < currentRow.size(); elementIdx++)
            {
                double sum = 0.0;
                int j = currentRow[elementIdx].first;
                for(int innerElementIdx = 1; innerElementIdx < currentRow.size(); innerElementIdx++)
                {
                    int i = currentRow[innerElementIdx].first;
                    if(output.isStored(j,i))
                    {
                        double temp = output.getValue(j,i) -
                                output.getValue(rowIdx,i) * output.getValue(rowIdx,j);
                        output.setValue(j,i,temp);
                    }
                    else
                    {
                        sum += output.getValue(rowIdx,i) * output.getValue(rowIdx,j);
                    }
                }
                output.setValue(j,j,output.getValue(j,j) - tuning * sum);
            }
        }
    }

    return output;
}

void LinearSolver::updateSubgrids()
{
    propagateMaterialGrid(m_mainMaterialGrid, m_materialSubgrids[0]);
    for(int i = 0; i < (m_materialSubgrids.size() - 1); i++)
    {
        propagateMaterialGrid(m_materialSubgrids[i], m_materialSubgrids[i+1]);
    }
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
                    coarseGrid.at(i,j) = FluidMaterial::EMPTY;
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

void LinearSolver::restrictGrid(const MaterialGrid &coarseMaterials, const Grid2d<double> &fineGrid, Grid2d<double> coarseGrid)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(coarseMaterials.data().size(),128);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::restrictGridThread,this,range,coarseMaterials,
                                 std::ref(fineGrid),std::ref(coarseGrid));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    ThreadPool::i()->wait();
}

void LinearSolver::restrictGridThread(const Range range, const MaterialGrid &coarseMaterials, const Grid2d<double> &fineGrid, Grid2d<double> coarseGrid)
{
    const std::vector<double> &fineGridData = fineGrid.data();
    const int gridDataSize = fineGridData.size();
    for(int idx = range.start; idx < range.end; idx++)
    {
        Index2d i2d = coarseMaterials.index2d(idx);
        std::array<int,16> stencilCellIdxs = getFineGridStencilIdxs(fineGrid,i2d.m_i,i2d.m_j);
        coarseGrid.at(i2d.m_i,i2d.m_j) = 0.f;
        if(!coarseMaterials.isFluid(i2d.m_i,i2d.m_j))
        {
            continue;
        }

        int weightIdx = 0;
        for(int idx : stencilCellIdxs)
        {
            if(idx < 0 || idx >= gridDataSize)
            {
                weightIdx++;
                continue;
            }
            coarseGrid.at(i2d.m_i,i2d.m_j) += fineGridData[idx] * m_restrictionWeights[weightIdx];
            weightIdx++;
        }
    }
}

void LinearSolver::prolongateGrid(const MaterialGrid &fineMaterials, const Grid2d<double> &coarseGrid, std::vector<double>& fineGridData)
{
    std::vector<Range> ranges = ThreadPool::i()->splitRange(fineMaterials.data().size(),128);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::prolongateGridThread,this,range,fineMaterials,
                                 std::ref(coarseGrid),std::ref(fineGridData));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    ThreadPool::i()->wait();
}

void LinearSolver::prolongateGridThread(const Range range, const MaterialGrid &fineMaterials, const Grid2d<double> &coarseGrid, std::vector<double> &fineGridData)
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
        std::array<int, 4> coarseNeghbors = getCoarseProlongIdxs(coarseGrid,i2d.m_i,i2d.m_j);
        std::array<float,4> coarseWeights = getProlongationWeights(i2d.m_i,i2d.m_j);
        for(int idx = 0; idx < 4; idx++)
        {
            if(idx < 0 || idx >= gridDataSize)
            {
                continue;
            }
            fineGridData[fineMaterials.linearIndex(i2d.m_i,i2d.m_j)] +=
                coarseGridData[coarseNeghbors[idx]] * coarseWeights[idx];
        }
    }
}

void LinearSolver::dampedJacobi(const MaterialGrid &materials, std::vector<double> &pressures, const std::vector<double> &rhs)
{
    std::vector<double> temp(pressures.size(), 0.);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(pressures.size(),128);
    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(m_dampedJacobiThread,this,range,std::ref(materials),
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
    const float tune = 2.f/3.f;
    for(int i = range.start; i < range.end; i++)
    {
        const auto weights = solver->getMultigridMatrixEntriesForCell(materials,i);
        int currIdx = weights[0].first;
        float result = 0.f;
        for(int wIdx = 1; wIdx < weights.size(); wIdx++)
        {
            if(weights[wIdx].first < 0)
            {
                continue;
            }

            result += weights[wIdx].second * pressures[weights[wIdx].first];
        }
        vout[currIdx] = ((rhs[currIdx]-result)/weights[0].second) * tune;
    }
}

void LinearSolver::vaddmul(const Range range, std::vector<double> &vin, const std::vector<double> &vadd, double weight)
{
    for(int i = range.start; i < range.end; i++)
    {
        vin[i] += vadd[i] * weight;
    }
}

void LinearSolver::vCycle(std::vector<double> &vout, const std::vector<double> &vin)
{
    vout.assign(vin.size(),0.f);
    const int subLevelCount = m_materialSubgrids.size();
    //Down stroke
    for(int level = -1; level < static_cast<int>(m_materialSubgrids.size() - 1); level++)
    {
        if(level != -1)
        {
            Grid2d<double> temp(m_materialSubgrids[level].sizeI(),m_materialSubgrids[level].sizeJ(),0.f);
            dampedJacobi(m_materialSubgrids[level],m_pressureGrids[level].data(),m_rhsGrids[level].data());
            multigridSubMatmul(m_materialSubgrids[level],
                               m_rhsGrids[level].data(),
                               m_pressureGrids[level].data(),
                               temp.data());
            restrictGrid(m_materialSubgrids[level+1],temp,m_rhsGrids[level+1]);
            m_pressureGrids[level+1].data().assign(m_pressureGrids[level+1].data().size(),0.f);
        }
        else
        {
            Grid2d<double> temp(m_mainMaterialGrid.sizeI(),m_mainMaterialGrid.sizeJ(),0.f);
            dampedJacobi(m_mainMaterialGrid,vout,vin);
            multigridSubMatmul(m_mainMaterialGrid,
                               vin,
                               vout,
                               temp.data());
            restrictGrid(m_materialSubgrids[0],temp,m_rhsGrids[0]);
            m_pressureGrids[0].data().assign(m_pressureGrids[0].data().size(),0.f);
        }
    }

    //Solve coarse
    for(int i = 0; i < 10; i++)
    {
        dampedJacobi(m_materialSubgrids[subLevelCount-1],
                     m_pressureGrids[subLevelCount-1].data(),
                     m_rhsGrids[subLevelCount-1].data());
    }


    //Up stroke
    for(int level = m_materialSubgrids.size() - 2; level >= -1; level--)
    {
        if(level != -1)
        {
            prolongateGrid(m_materialSubgrids[level],m_pressureGrids[level+1],m_pressureGrids[level].data());
            dampedJacobi(m_materialSubgrids[level],m_pressureGrids[level].data(),m_rhsGrids[level].data());
        }
        else
        {
            prolongateGrid(m_mainMaterialGrid,m_pressureGrids[0],vout);
            dampedJacobi(m_mainMaterialGrid,vout,vin);
        }
    }
}

void LinearSolver::multigridMatmul(const MaterialGrid &materials, const std::vector<double> &vin, std::vector<double> &vout)
{
    for(int idx = 0; idx < vin.size(); idx++)
    {
        const auto weights = getMultigridMatrixEntriesForCell(materials,idx);
        float result = 0.f;
        for(const auto& w : weights)
        {
            if(w.first < 0)
            {
                continue;
            }

            result += w.second * vin[w.first];
        }
        vout[idx] = result;
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
    for(int idx = range.start; idx < range.end; idx++)
    {
        const auto weights = getMultigridMatrixEntriesForCell(materials,idx);
        float result = 0.f;
        for(const auto& w : weights)
        {
            if(w.first < 0)
            {
                continue;
            }

            result += w.second * vmul[w.first];
        }
        vout[idx] = vsub[idx] - result;
    }
}

std::array<std::pair<int, float>, 5> LinearSolver::getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int i, int j)
{
    std::array<std::pair<int,float>,5> output{std::pair<int,float>()};
    std::array<int,4> neighbors = materials.immidiateNeighbors(i,j);
    const std::vector<FluidMaterial> &materialsData = materials.data();

    int neighborCount = 0;
    int outputIdx = 1;
    for(int idx : neighbors)
    {
        if(idx < 0 || idx > materialsData.size())
        {
            output[outputIdx] = std::make_pair(-1,0);
            outputIdx++;
            continue;
        }
        const float weight = !solidTest(materialsData[idx]);
        output[outputIdx] = std::make_pair(idx,weight);
        neighborCount++;
        outputIdx++;
    }
    output[0] = std::make_pair(materials.linearIndex(i,j),-neighborCount);

    return output;
}

std::array<std::pair<int, float>, 5> LinearSolver::getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int linearIdx)
{
    std::array<std::pair<int,float>,5> output{std::pair<int,float>()};
    std::array<int,4> neighbors = materials.immidiateNeighbors(linearIdx);
    const std::vector<FluidMaterial> &materialsData = materials.data();
    const auto dataSize = materialsData.size();

    int neighborCount = 0;
    int outputIdx = 1;
    for(int idx : neighbors)
    {
        if(idx < 0 || idx >= dataSize)
        {
            output[outputIdx] = std::make_pair(-1,0);
            outputIdx++;
            continue;
        }
        const float weight = solidTest(materialsData[idx]) ? 0.f : 1.f;
        output[outputIdx] = std::make_pair(idx,weight);
        neighborCount++;
        outputIdx++;
    }
    output[0] = std::make_pair(linearIdx,-neighborCount);

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
