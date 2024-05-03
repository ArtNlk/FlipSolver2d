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
#include "pressuredata.h"
#include "vmath.h"
#include "linearsolver_sse42.h"

#include "logger.h"

LinearSolver::LinearSolver()
{

}

bool LinearSolver::solve(const IndexedPressureParameters &matrixIn,
                         std::vector<double> &result,
                         const std::vector<double> &vec,
                         int iterLimit,
                         double tol)
{
   result.assign(result.size(),0);
   if (VOps::i().isZero(vec))
   {
       std::cout << "Solver skipped zeros vector" << '\n';
       return true;
   }
   //DynamicUpperTriangularSparseMatrix precond = calcPrecond(matrixIn);
   //UpperTriangularMatrix matrix(matrixIn);

   //debug() << "mat=" << matrix;
   std::vector<double> residual(vec);
   std::vector<double> aux = residual;
   //applyICPrecond(precond,residual,aux);
   //applyIPPrecond(&matrix,&residual,&aux);
   //precond->apply(residual,aux,data);
   std::vector<double> search = aux;
   double sigma = VOps::i().dot(aux,residual);
   double err = 0.0;
   for (int i = 0; i < iterLimit; i++)
   {
       matrixIn.multiply(search, aux);
       double alpha = sigma/(VOps::i().dot(aux,search) + 1e-8);
       VOps::i().addMul(result,result,search,alpha);
       VOps::i().subMul(residual,residual,aux,alpha);
       err = VOps::i().maxAbs(residual);
       if(i % 5 == 0)
       {
           std::cout << "Solver: " << i << " : " << err << "\n";
           debug() << "Solver: " << i << " : " << err;
       }
       if (err <= tol)
       {
           debug() << "[SOLVER] Solver done, iter = " << i << " err = " << err;
           std::cout << "Solver done, iter = " << i << " err = " << err << '\n';
           return true;
       }
       aux = residual;
       //applyICPrecond(precond,residual,aux);
       //applyIPPrecond(&matrix,&residual,&aux);
       double newSigma = VOps::i().dot(aux,residual);
       double beta = newSigma/(sigma);
       VOps::i().addMul(search,aux,search,beta);
       sigma = newSigma;
   }

   debug() << "Solver iter exhaustion, err = " << err;
   std::cout << "Solver iter exhaustion, err = " << err << '\n';
    return false;
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
