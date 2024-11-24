#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <vector>

#include "pressuredata.h"
#include "PressureIPPCoeficients.h"

class LinearSolver
{
public:
    LinearSolver(MaterialGrid &materialGrid, int maxMultigridDepth);
    using SparseMatRowElements = std::array<std::pair<int,double>,5>;

    int solve(const IndexedPressureParameters &matrixIn,
               const IndexedIPPCoefficients &precond,
               std::vector<double> &result,
               const std::vector<double> &vec,
               int iterLimit = 20,
               double tol = 1e-6);

    friend class LinearSolver_sse42;

protected:
    void applyMGPrecond(std::vector<double> const &in, std::vector<double> &out);

    void updateSubgrids();

    void propagateMaterialGrid(const MaterialGrid& fineGrid, MaterialGrid &coarseGrid);

    void restrictGrid(const MaterialGrid& coarseMaterials, const Grid2d<double> &fineGrid,
                      Grid2d<double> coarseGrid);

    void restrictGridThread(const Range range, const MaterialGrid& coarseMaterials,
                            const Grid2d<double> &fineGrid, Grid2d<double> coarseGrid);
    void prolongateGrid(const MaterialGrid& fineMaterials, const Grid2d<double> &coarseGrid,
                        std::vector<double> &fineGridData);

    void prolongateGridThread(const Range range, const MaterialGrid& fineMaterials,
                              const Grid2d<double> &coarseGrid, std::vector<double> &fineGridData);
    void premaskPressures(const MaterialGrid& materials, std::vector<double>& pressures);

    static void premaskPressuresThread(const Range range, const MaterialGrid& materials, std::vector<double>& pressures);

    void dampedJacobi(const MaterialGrid& materials, std::vector<double> &pressures, const std::vector<double> &rhs);

    static void dampedJacobiThread(LinearSolver* solver, const Range range, const MaterialGrid& materials, std::vector<double> &vout,
                                   const std::vector<double> &pressures, const std::vector<double> &rhs);

    void vCycle(std::vector<double> &vout, const std::vector<double> &vin);
    void multigridSubMatmul(const MaterialGrid& materials, const std::vector<double>& vsub,
                            const std::vector<double>& vmul, std::vector<double>& vout);
    void multigridSubMatmulThread(const Range range, const MaterialGrid& materials, const std::vector<double>& vsub,
                                  const std::vector<double>& vmul, std::vector<double>& vout);

    std::pair<std::array<int,4>,std::array<double,4>> getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int i, int j);
    std::pair<std::array<int,4>,std::array<double,4>> getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int linearIdx);
    std::array<int,16> getFineGridStencilIdxs(const LinearIndexable2d &fineGridIndexer,int iCoarse, int jCoarse);
    std::array<int,4> getFineGridChildIdxs(const LinearIndexable2d &fineGridIndexer,int iCoarse, int jCoarse);
    std::array<float,4> getProlongationWeights(int finei, int finej);
    std::array<int,4> getCoarseProlongIdxs(const LinearIndexable2d &coarseGridIndexer, int finei, int finej);
    std::array<float,16> m_restrictionWeights;
    std::array<float,16> m_prolongationWeights;

    MaterialGrid& m_mainMaterialGrid;
    std::vector<MaterialGrid> m_materialSubgrids;
    std::vector<Grid2d<double>> m_pressureGrids;
    std::vector<Grid2d<double>> m_rhsGrids;

};
#endif // PCGSOLVER_H
