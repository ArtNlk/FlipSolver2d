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
#include "vmath.h"

#include "logger.h"

const double LinearSolver::m_tol = 1.0e-6;

LinearSolver::LinearSolver(MaterialGrid &materialGrid, int maxMultigridDepth) :
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
        m_materialSubgrids.emplace_back(sizeI,sizeJ,FluidMaterial::SOLID);
        m_pressureGrids.emplace_back(sizeI,sizeJ,0.f,OOBStrategy::OOB_EXTEND);
        m_rhsGrids.emplace_back(sizeI,sizeJ,0.f,OOBStrategy::OOB_EXTEND);
    }

    const std::array<float, 4> weights = {1.f,3.f,3.f,1.f};

    int idx = 0;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            m_restrictionWeights[idx] = (weights[i] * weights[j]) / 64.f;
            m_prolongationWeights[idx] = weights[i] * weights[j];
        }
    }
}

bool LinearSolver::solve(const DynamicUpperTriangularSparseMatrix &matrixIn, std::vector<double> &result, const std::vector<double> &vec, int iterLimit)
{
    result.assign(result.size(),0);
    if (vsimmath::isZero(vec))
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
    double sigma = vsimmath::dot(aux,residual);
    double err = 0.0;
    for (int i = 0; i < iterLimit; i++)
    {
        aux = matrix * search;
        double alpha = sigma/(vsimmath::dot(aux,search));
        vsimmath::addMul(result,result,search,alpha);
        vsimmath::subMul(residual,residual,aux,alpha);
        err = vsimmath::maxAbs(residual);
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
        double newSigma = vsimmath::dot(aux,residual);
        double beta = newSigma/(sigma);
        vsimmath::addMul(search,aux,search,beta);
        sigma = newSigma;
    }

    debug() << "Solver iter exhaustion, err = " << err;
    std::cout << "Solver iter exhaustion, err = " << err << '\n';
    return false;
}

bool LinearSolver::mfcgSolve(MatElementProvider elementProvider, std::vector<double> &result, const std::vector<double> &vec, int iterLimit)
{
    result.assign(result.size(),0);
    if (vsimmath::isZero(vec))
    {
        return true;
    }

    //vec - b
    std::vector<double> residual = vec; //r
    std::vector<double> aux = vec; //q
    //applyMGPrecond(residual,aux);
    std::vector<double> search = aux; //d
    double sigma = vsimmath::dot(aux,residual);
    double err = 0.0;
    for (int i = 0; i < iterLimit; i++)
    {
        nomatVMul(elementProvider,search,aux);
        double alpha = sigma/(vsimmath::dot(aux,search));
        vsimmath::addMul(result,result,search,alpha);
        vsimmath::subMul(residual,residual,aux,alpha);
        err = vsimmath::maxAbs(residual);
        if (err <= m_tol)
        {
            //debug() << "[SOLVER] Solver done, iter = " << i << " err = " << err;
            std::cout << "Solver done, iter = " << i << " err = " << err << '\n';
            return true;
        }
        aux = residual;
        applyMGPrecond(residual,aux);
        double newSigma = vsimmath::dot(aux,residual);
        double beta = newSigma/(sigma);
        vsimmath::addMul(search,aux,search,beta);
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
    const std::vector<double> &fineGridData = fineGrid.data();
    const int gridDataSize = fineGridData.size();
    for(int i = 0; i < coarseGrid.sizeI(); i++)
    {
        for(int j = 0; j < coarseGrid.sizeJ(); j++)
        {
            std::array<int,16> stencilCellIdxs = getFineGridStencilIdxs(fineGrid,i,j);
            coarseGrid.at(i,j) = 0.f;
            if(!coarseMaterials.isFluid(i,j))
            {
                continue;
            }

            for(int idx : stencilCellIdxs)
            {
                if(idx < 0 || idx >= gridDataSize)
                {
                    continue;
                }
                coarseGrid.at(i,j) += fineGridData[idx] * m_restrictionWeights[idx];
            }
        }
    }
}

void LinearSolver::prolongateGrid(const MaterialGrid &fineMaterials, const Grid2d<double> &coarseGrid, std::vector<double>& fineGridData)
{
    const std::vector<double> &coarseGridData = coarseGrid.data();
    const int gridDataSize = coarseGridData.size();
    for(int i = 0; i < fineMaterials.sizeI(); i++)
    {
        for(int j = 0; j < fineMaterials.sizeJ(); j++)
        {
            if(!fineMaterials.isFluid(i,j))
            {
                continue;
            }
            std::array<int, 4> coarseNeghbors = getCoarseProlongIdxs(coarseGrid,i,j);
            std::array<float,4> coarseWeights = getProlongationWeights(i,j);
            for(int idx = 0; idx < 4; idx++)
            {
                if(idx < 0 || idx >= gridDataSize)
                {
                    continue;
                }
                fineGridData[fineMaterials.linearIndex(i,j)] +=
                    coarseGridData[coarseNeghbors[idx]] * coarseWeights[idx];
            }
        }
    }
}

void LinearSolver::dampedJacobi(const MaterialGrid &materials, std::vector<double> &pressures, const std::vector<double> &rhs)
{
    std::vector<double> temp(pressures.size(),0.f);
    const float tune = 2.f/3.f;
    for(int i = 0; i < pressures.size(); i++)
    {
        const auto weights = getMultigridMatrixEntriesForCell(materials,i);
        int currIdx = weights[0].first;
        float result = 0.f;
        for(const auto& w : weights)
        {
            if(w.first < 0)
            {
                continue;
            }

            result += w.second * pressures[w.first];
        }
        temp[currIdx] = (rhs[currIdx]-result)/weights[0].second;
    }

    for(int i = 0; i < pressures.size(); i++)
    {
        pressures[i] = temp[i] * tune;
    }
}

void LinearSolver::vCycle(std::vector<double> &vout, const std::vector<double> &vin)
{
    vout.assign(vin.size(),0.f);
    const int subLevelCount = m_materialSubgrids.size();
    //Down stroke
    for(int level = -1; level < m_materialSubgrids.size() - 1; level++)
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
        if(level > -1)
        {
            prolongateGrid(m_materialSubgrids[level],m_pressureGrids[level+1],m_pressureGrids[level].data());
            dampedJacobi(m_materialSubgrids[level],m_pressureGrids[level].data(),m_rhsGrids[level].data());
        }
        else
        {
            prolongateGrid(m_mainMaterialGrid,m_pressureGrids[level+1],vout);
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
    for(int idx = 0; idx < vout.size(); idx++)
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
