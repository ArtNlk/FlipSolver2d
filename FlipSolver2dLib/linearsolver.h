#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <vector>

#include "PressureWeights.h"
#include "InversePoissonPreconditioner.h"

#include "logger.h"
#include "vmath.h"

class LinearSolver
{
public:
    LinearSolver() = default;

    template<class T, class U>
    int solve(const MatrixWeights<T> &matrixIn,
               const MatrixWeights<U> &precond,
               std::vector<double> &result,
               const std::vector<double> &vec,
               int iterLimit = 20,
               double tol = 1e-6)
    {
        result.assign(result.size(), 0);
        if (VOps::i().isZero(vec)) {
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
        double sigma = VOps::i().dot(aux, residual);
        double err = 0.0;
        for (int i = 1; i <= iterLimit; i++) {
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
                debug() << "[SOLVER] Solver done, iter = " << i << " err = " << err;
                std::cout << "Solver done, iter = " << i << " err = " << err << '\n';

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

        debug() << "Solver iter exhaustion, err = " << err;
        std::cout << "Solver iter exhaustion, err = " << err << '\n';
        return iterLimit;
    }

    friend class LinearSolver_sse42;
};
#endif // PCGSOLVER_H
