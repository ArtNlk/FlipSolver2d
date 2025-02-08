#include "linearsolver.h"

#include <algorithm>

#include <array>
#include <cmath>
#include <future>
#include <mutex>
#include <sstream>
#include <utility>
#include <vector>

#include "PressureIPPCoeficients.h"
#include "dynamicmatrix.h"
#include "grid2d.h"
#include "materialgrid.h"
#include "pressuredata.h"
#include "vmath.h"
#include "linearsolver_sse42.h"

#include <fstream>

#include "logger.h"

int LinearSolver::solve(const IndexedPressureParameters &matrixIn,
                        const IndexedIPPCoefficients &precond,
                        std::vector<double> &result,
                        const std::vector<double> &vec,
                        int iterLimit,
                        double tol)
{
    result.assign(result.size(), 0);
    if (VOps::i().isZero(vec)) {
        return 0;
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
    double sigma = VOps::i().dot(aux, residual);
    double err = 0.0;
    for (int i = 0; i < iterLimit; i++) {
        matrixIn.multiply(search, aux);
        double alpha = sigma / (VOps::i().dot(aux, search) + 1e-8);
        VOps::i().addMul(result, result, search, alpha);
        VOps::i().subMul(residual, residual, aux, alpha);
        err = VOps::i().maxAbs(residual);
        // if(i % 5 == 0)
        // {
        //     std::cout << "Solver: " << i << " : " << err << "\n";
        //     debug() << "Solver: " << i << " : " << err;
        // }
        if (err <= tol) {
            return i;
        }
        //aux = residual;
        precond.multiply(residual, aux);
        //applyICPrecond(precond,residual,aux);
        //applyIPPrecond(&matrix,&residual,&aux);
        double newSigma = VOps::i().dot(aux, residual);
        double beta = newSigma / (sigma);
        VOps::i().addMul(search, aux, search, beta);
        sigma = newSigma;
    }

    return iterLimit;
}
