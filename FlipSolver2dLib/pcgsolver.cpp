#include "pcgsolver.h"

#include <algorithm>

#include <sstream>
#include <vector>

#include "dynamicuppertriangularsparsematrix.h"
#include "vmath.h"
#include "cmath"
#include "simsettings.h"
#include "logger.h"

const double PCGSolver::m_tol = 1.0e-14;

bool PCGSolver::solve(const DynamicUpperTriangularSparseMatrix &matrixIn, std::vector<double> &result, const std::vector<double> &vec, int iterLimit)
{
    result.assign(result.size(),0);
    if (vmath::isZero(vec))
    {
        return true;
    }
    DynamicUpperTriangularSparseMatrix precond = calcPrecond(matrixIn);
    UpperTriangularMatrix matrix(matrixIn);
    std::vector<double> residual = vec;
    std::vector<double> aux = vec;
    std::vector<double> aux_before = aux;
    applyICPrecond(precond,residual,aux);
    std::vector<double> search = aux;
    double sigma = vmath::dot(aux,residual);
    double err = 0.0;
    for (int i = 0; i < iterLimit; i++)
    {
        aux = matrix * search;
        double alpha = sigma/(vmath::dot(aux,search));
        vmath::addMul(result,result,search,alpha);
        vmath::subMul(residual,residual,aux,alpha);
        err = vmath::maxAbs(residual);
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
        double newSigma = vmath::dot(aux,residual);
        double beta = newSigma/(sigma);
        vmath::addMul(search,aux,search,beta);
        sigma = newSigma;
    }

    debug() << "Solver iter exhaustion, err = " << err;
    std::cout << "Solver iter exhaustion, err = " << err << '\n';
    return false;
}

void PCGSolver::applyICPrecond(const DynamicUpperTriangularSparseMatrix &precond, const std::vector<double> &in, std::vector<double> &out)
{
    out = std::vector(in);

    std::vector<SparseRow> &rows = const_cast<DynamicUpperTriangularSparseMatrix&>(precond).data();

    for(int i = 0; i < out.size(); i++)
    {
        auto& currentRow = rows[i];
        if(precond.getValue(i,i) != 0.0)
        {
            out[i] /= precond.getValue(i,i);

            for(int elementIdx = 1; elementIdx < currentRow.size(); elementIdx++)
            {
                int j = currentRow[elementIdx].first;
                out[j] = out[j] - precond.getValue(i,j) * out[j];
            }
        }
    }

    for(int i = out.size() - 1; i >= 0; i--)
    {
        auto& currentRow = rows[i];
        if(precond.getValue(i,i) != 0.0)
        {
            for(int elementIdx = 1; elementIdx < currentRow.size(); elementIdx++)
            {
                int j = currentRow[elementIdx].first;
                out[i] = out[i] - precond.getValue(i,j) * out[j];
            }
            out[i] *= precond.getValue(i,i);
        }
    }
}

DynamicUpperTriangularSparseMatrix PCGSolver::calcPrecond(const DynamicUpperTriangularSparseMatrix &matrix)
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
            double diag = currentRow[0].second;
            currentRow[0].second = sqrt(diag);

            for(int elementIdx = 1; elementIdx < currentRow.size(); elementIdx++)
            {
                currentRow[elementIdx].second = currentRow[elementIdx].second / diag;
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
