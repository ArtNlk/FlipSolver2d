#include "linearsolver.h"

#include <algorithm>

#include <array>
#include <cmath>
#include <future>
#include <mutex>
#include <sstream>
#include <utility>
#include <vector>

#include "dynamicmatrix.h"
#include "grid2d.h"
#include "materialgrid.h"
#include "vmath.h"
#include "linearsolver_sse42.h"

#include "logger.h"

LinearSolver::LinearSolver(MaterialGrid &materialGrid, int maxMultigridDepth) :
    m_restrictionWeights({0}),
    m_prolongationWeights({0}),
    m_mainMaterialGrid(materialGrid),
    m_premaskPressuresThread(&LinearSolver::premaskPressuresThread)
{
#ifdef FLUID_SSE
        m_premaskPressuresThread = &LinearSolver_sse42::premaskPressuresThread;
#endif

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

bool LinearSolver::solve(const StaticMatrix &matrixIn, std::vector<double> &result, const std::vector<double> &vec, std::shared_ptr<IPreconditioner> precond, IPreconditioner::PreconditionerData* data, int iterLimit, double tol)
{
//    result.assign(result.size(),0);
//    if(precond == nullptr)
//    {
//        return false;
//    }
//    if (VOps::i().isZero(vec))
//    {
//        std::cout << "Solver skipped zeros vector" << '\n';
//        return true;
//    }
//    //DynamicUpperTriangularSparseMatrix precond = calcPrecond(matrixIn);
//    //UpperTriangularMatrix matrix(matrixIn);

//    //debug() << "mat=" << matrix;
//    Eigen::VectorXd residual(vec.data());
//    Eigen::VectorXd aux = residual;
//    //applyICPrecond(precond,residual,aux);
//    //applyIPPrecond(&matrix,&residual,&aux);
//    //precond->apply(residual,aux,data);
//    Eigen::VectorXd search = aux;
//    double sigma = VOps::i().dot(aux,residual);
//    double err = 0.0;
//    for (int i = 0; i < iterLimit; i++)
//    {
//        aux = matrixIn * search;
//        double alpha = sigma/(VOps::i().dot(aux,search) + 1e-8);
//        VOps::i().addMul(result,result,search,alpha);
//        VOps::i().subMul(residual,residual,aux,alpha);
//        err = VOps::i().maxAbs(residual);
////        if(i % 5 == 0)
////        {
////            std::cout << "Solver: " << i << " : " << err << "\n";
////            debug() << "Solver: " << i << " : " << err;
////        }
//        if (err <= tol)
//        {
//            debug() << "[SOLVER] Solver done, iter = " << i << " err = " << err;
//            std::cout << "Solver done, iter = " << i << " err = " << err << '\n';
//            return true;
//        }
//        //aux = residual;
//        //applyICPrecond(precond,residual,aux);
//        //applyIPPrecond(&matrix,&residual,&aux);
//        precond->apply(residual,aux,data);
//        double newSigma = VOps::i().dot(aux,residual);
//        double beta = newSigma/(sigma);
//        VOps::i().addMul(search,aux,search,beta);
//        sigma = newSigma;
//    }

//    debug() << "Solver iter exhaustion, err = " << err;
//    std::cout << "Solver iter exhaustion, err = " << err << '\n';
    return false;
}

bool LinearSolver::mfcgSolve(MatElementProvider elementProvider, std::vector<double> &result, const std::vector<double> &vec, std::shared_ptr<IPreconditioner> precond, int iterLimit, double tol)
{
//    result.assign(result.size(),0);
//    if(precond == nullptr)
//    {
//        return false;
//    }
//    if (VOps::i().isZero(vec))
//    {
//        return true;
//    }
//    updateSubgrids();

//    //vec - b
//    std::vector<double> residual = vec; //r
//    std::vector<double> aux = vec; //z
//    //applyMGPrecond(residual,aux);
//    //applyMfIPPrecond(elementProvider,residual,aux);
//    precond->apply(residual,aux);
//    std::vector<double> search = aux; //p
//    double sigma = VOps::i().dot(aux,residual);
//    double err = 0.0;
//    for (int i = 0; i < iterLimit; i++)
//    {
//        nomatVMul(elementProvider,search,aux);
//        double alpha = sigma/(VOps::i().dot(aux,search));
//        VOps::i().subMul(residual,residual,aux,alpha);
//        err = VOps::i().maxAbs(residual);
//        if (err <= tol)
//        {
//            //debug() << "[SOLVER] Solver done, iter = " << i << " err = " << err;
//            std::cout << "MFSolver done, iter = " << i << " err = " << err << '\n';
//            return true;
//        }
//        //aux = residual;
//        //applyMGPrecond(residual,aux);
//        //applyMfIPPrecond(elementProvider,residual,aux);
//        precond->apply(residual,aux);
//        double newSigma = VOps::i().dot(aux,residual);
//        //std::cout << "New sigma:" << newSigma << std::endl;
//        double beta = newSigma/(sigma);
//        VOps::i().addMul(result,result,search,alpha);
//        VOps::i().addMul(search,aux,search,beta);
//        sigma = newSigma;
//    }

//    //debug() << "Solver iter exhaustion, err = " << err;
//    //std::cout << "Solver iter exhaustion, err = " << err << '\n';
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

void LinearSolver::applyICPrecond(const DynamicMatrix &precond, const std::vector<double> &in, std::vector<double> &out)
{
    out = in;
    
    std::vector<SparseRow> &rows = const_cast<DynamicMatrix&>(precond).data();

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

DynamicMatrix LinearSolver::calcPrecond(const DynamicMatrix &matrix)
{
    float tuning = 0.97f;
    float safety = 0.25f;
    int n = matrix.size();
    DynamicMatrix output(n);
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

void LinearSolver::applyMfIPPrecond(MatElementProvider p, const std::vector<double> &in, std::vector<double> &out)
{
//    out = in;
//    return;
    std::vector<double> temp(in.size(),0.0);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(in.size(),128);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::firstStepMfIPPMatmulThread,range,p,
                                 std::ref(in),std::ref(temp));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    ThreadPool::i()->wait();

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&LinearSolver::secondStepMfIPPMatmulThread,range,p,
                                 std::ref(temp),std::ref(out));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    ThreadPool::i()->wait();

}

void LinearSolver::firstStepMfIPPMatmulThread(Range r, MatElementProvider p, const std::vector<double> &in, std::vector<double> &out)
{
    for(int i = r.start; i < r.end; i++)
    {
        SparseMatRowElements matElements = p(i);
        if(std::abs(matElements[4].second) < 0.0001)
        {
            out[i] = in[i];
            continue;
        }
        out[i] = in[i] - (in[matElements[1].first] * matElements[1].second
                           + in[matElements[3].first] * matElements[3].second) / matElements[4].second;
    }
}

void LinearSolver::secondStepMfIPPMatmulThread(Range r, MatElementProvider p, const std::vector<double> &in, std::vector<double> &out)
{
    for(int i = r.start; i < r.end; i++)
    {
        SparseMatRowElements matElements = p(i);
        if(std::abs(matElements[4].second) < 0.0001)
        {
            out[i] = in[i];
            continue;
        }
        out[i] = in[i] - (in[matElements[0].first] * matElements[0].second
                            + in[matElements[2].first] * matElements[2].second) / matElements[4].second;;
    }
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
        std::array<int,16> stencilCellIdxs = getFineGridStencilIdxs(fineGrid,i2d.i,i2d.j);
        coarseGrid.at(i2d.i,i2d.j) = 0.f;
        if(!coarseMaterials.isFluid(i2d.i,i2d.j))
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
            coarseGrid.at(i2d.i,i2d.j) += fineGridData[idx] * m_restrictionWeights[weightIdx];
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
        std::array<int, 4> coarseNeghbors = getCoarseProlongIdxs(coarseGrid,i2d.i,i2d.j);
        std::array<float,4> coarseWeights = getProlongationWeights(i2d.i,i2d.j);
        for(int idx = 0; idx < 4; idx++)
        {
            if(idx < 0 || idx >= gridDataSize)
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
        ThreadPool::i()->enqueue(m_premaskPressuresThread,range,
                                 std::ref(materials),std::ref(pressures));
    }
    ThreadPool::i()->wait();
}

void LinearSolver::premaskPressuresThread(const Range range, const MaterialGrid &materials, std::vector<double> &pressures)
{
    for(int i = range.start; i < range.end; i++)
    {
        pressures[i] *= !solidTest(materials.data()[i]);
    }
}

void LinearSolver::dampedJacobi(const MaterialGrid &materials, std::vector<double> &pressures, const std::vector<double> &rhs)
{
    std::vector<double> temp(pressures.size(), 0.);

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
    const double tune = 2.0/3.0;
    const size_t dataSize = materials.data().size();
    for(int i = range.start; i < range.end; i++)
    {
        //const auto weights = solver->getMultigridMatrixEntriesForCell(materials,i);
        const auto neighbors = materials.immidiateNeighbors(i);
        double result = 0.0;
        for(int idx : neighbors)
        {
            if(idx < 0 || idx >= dataSize)
            {
                continue;
            }
            result += pressures[idx];
        }
        vout[i] = ((rhs[i]-result)/-4.0) * tune;
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
            premaskPressures(m_materialSubgrids[level], m_pressureGrids[level].data());
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
            premaskPressures(m_mainMaterialGrid, vout);
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
        premaskPressures(m_materialSubgrids[subLevelCount-1], m_pressureGrids[subLevelCount-1].data());
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
            premaskPressures(m_materialSubgrids[level], m_pressureGrids[level].data());
            dampedJacobi(m_materialSubgrids[level],m_pressureGrids[level].data(),m_rhsGrids[level].data());
        }
        else
        {
            prolongateGrid(m_mainMaterialGrid,m_pressureGrids[0],vout);
            premaskPressures(m_mainMaterialGrid, vout);
            dampedJacobi(m_mainMaterialGrid,vout,vin);
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
        for(int neighborIdx : neighbors)
        {
            if(neighborIdx < 0 || neighborIdx >= dataSize)
            {
                continue;
            }
            result += vmul[neighborIdx];
        }
        result += -4.0 * vmul[idx];
        vout[idx] = vsub[idx] - result;
    }
}

std::pair<std::array<int,4>,std::array<double,4>> LinearSolver::getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int i, int j)
{
    return getMultigridMatrixEntriesForCell(materials, materials.linearIndex(i,j));
}

std::pair<std::array<int,4>,std::array<double,4>> LinearSolver::getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int linearIdx)
{
    std::pair<std::array<int,4>,std::array<double,4>> output;
    std::array<int,4> neighbors = materials.immidiateNeighbors(linearIdx);
    const std::vector<FluidMaterial> &materialsData = materials.data();
    const auto dataSize = materialsData.size();

    int outputIdx = 0;
    for(int idx : neighbors)
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

void StubPreconditioner::apply(const std::vector<double> &in, std::vector<double> &out, PreconditionerData *data)
{
    out.assign(in.begin(),in.end());
}

IPPreconditioner::IPPreconditioner(LinearIndexable2d indexer) :
    m_indexer(indexer)
{

}

void IPPreconditioner::apply(const std::vector<double> &in, std::vector<double> &out, PreconditionerData *data)
{
    if(data == nullptr)
    {
        std::cout << "IPP data null!" << std::endl;
        return;
    }

    IPPreconditionerData* ippData = dynamic_cast<IPPreconditionerData*>(data);
    std::vector<double> temp(in.size(),0.0);

    std::vector<Range> ranges = ThreadPool::i()->splitRange(in.size(),128);

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&IPPreconditioner::firstStepIPPMatmulThread,this,range,
                                 std::ref(in),std::ref(temp),std::ref(ippData->m));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    ThreadPool::i()->wait();

    for(Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&IPPreconditioner::secondStepIPPMatmulThread,this,range,
                                 std::ref(temp),std::ref(out),std::ref(ippData->m));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    ThreadPool::i()->wait();
}

void IPPreconditioner::firstStepIPPMatmulThread(Range r, const std::vector<double> &in, std::vector<double> &out, StaticMatrix &m)
{
    for(int i = r.start; i < r.end; i++)
    {
        auto neighbors = m_indexer.immidiateNeighbors(i);
        double diag = m.getValue(i,i);
        if(std::abs(diag) < 0.0001)
        {
            out[i] = in[i];
            continue;
        }
        double v1 = in[neighbors[1]];
        double v2 = in[neighbors[3]];
        out[i] = in[i] - (v1 * m.getValue(i,neighbors[1])
                                          + v2 * m.getValue(i,neighbors[3])) / diag;
    }
}

void IPPreconditioner::secondStepIPPMatmulThread(Range r, const std::vector<double> &in, std::vector<double> &out, StaticMatrix &m)
{
    for(int i = r.start; i < r.end; i++)
    {
        auto neighbors = m_indexer.immidiateNeighbors(i);
        double diag = m.getValue(i,i);
        if(std::abs(diag) < 0.0001)
        {
            out[i] = in[i];
            continue;
        }
        out[i] = in[i] - (in[neighbors[0]] * m.getValue(i,neighbors[0])
                                          + in[neighbors[2]] * m.getValue(i,neighbors[2])) / diag;
    }
}
