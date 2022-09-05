#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <unordered_map>

#include "uppertriangularmatrix.h"
#include "fluidgrid.h"

class PCGSolver
{
public:
    PCGSolver();

    bool solve(const UpperTriangularMatrix &matrix, MACFluidGrid &grid, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 20);

protected:
    void applyICPrecond(const UpperTriangularMatrix &matrix, std::vector<double> &in, std::vector<double> &out, MACFluidGrid &grid);
    void calcPrecond(const UpperTriangularMatrix &matrix, MACFluidGrid &grid);
    double precond(const UpperTriangularMatrix &matrix, int i, int j, MACFluidGrid &grid);
    void reset();
    std::vector<double> m_residual;
    std::vector<double> m_aux;
    std::vector<double> m_search;
    std::unordered_map<int,double> m_precondCache;
    static const double m_tol;
};

#endif // PCGSOLVER_H
