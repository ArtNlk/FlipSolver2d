#include "linearsolver.h"

#include <algorithm>

#include <cmath>
#include <mutex>
#include <sstream>
#include <vector>

#include "dynamicuppertriangularsparsematrix.h"
#include "vmath.h"

#include "logger.h"

const double LinearSolver::m_tol = 1.0e-14;

LinearSolver::LinearSolver(ThreadPool &pool) :
    m_poolRef(pool)
{

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
            //std::cout << "Solver done, iter = " << i << " err = " << err << '\n';
            return true;
        }
        aux = residual;
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
    std::vector<Range> ranges = m_poolRef.splitRange(vin.size());

    for(Range range : ranges)
    {
        m_poolRef.enqueue(&LinearSolver::nomatVMulThread,this,range,elementProvider,
                          std::ref(vin),std::ref(vout));
        //nomatVMulThread(range,elementProvider,vin,vout);
    }
    m_poolRef.wait();

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
