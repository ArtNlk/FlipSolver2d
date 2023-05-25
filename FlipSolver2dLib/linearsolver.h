#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <functional>
#include <unordered_map>
#include <vector>

#include "linearindexable2d.h"
#include "materialgrid.h"
#include "simddispatcher.h"
#include "threadpool.h"
#include "uppertriangularmatrix.h"

class LinearSolver : SimdDispatcher
{
public:
    LinearSolver(MaterialGrid& materialGrid, int maxMultigridDepth);
    using SparseMatRowElements = std::array<std::pair<int,double>,5>;
    using MatElementProvider = std::function<SparseMatRowElements(int)>;

    bool solve(const DynamicUpperTriangularSparseMatrix &matrixIn, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 20);
    bool mfcgSolve(MatElementProvider elementProvider, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 20);

    friend class LinearSolver_sse42;

protected:
    void nomatVMul(MatElementProvider elementProvider, const std::vector<double> &vin, std::vector<double> &vout);
    void nomatVMulThread(Range range, MatElementProvider elementProvider,
                         const std::vector<double> &vin, std::vector<double> &vout);
    void applyMGPrecond(std::vector<double> const &in, std::vector<double> &out);
    void applyICPrecond(const DynamicUpperTriangularSparseMatrix &precond, std::vector<double> const &in, std::vector<double> &out);
    DynamicUpperTriangularSparseMatrix calcPrecond(const DynamicUpperTriangularSparseMatrix &matrix);

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

    void dampedJacobi(const MaterialGrid& materials, std::vector<double> &pressures, const std::vector<double> &rhs);

    static void dampedJacobiThread(LinearSolver* solver, const Range range, const MaterialGrid& materials, std::vector<double> &vout,
                            const std::vector<double> &pressures, const std::vector<double> &rhs);

    void vaddmul(const Range range, std::vector<double> &vin, const std::vector<double> &vadd, double weight);

    void vCycle(std::vector<double> &vout, const std::vector<double> &vin);

    void multigridMatmul(const MaterialGrid& materials, const std::vector<double>& vin, std::vector<double>& vout);

    void multigridSubMatmul(const MaterialGrid& materials, const std::vector<double>& vsub,
                            const std::vector<double>& vmul, std::vector<double>& vout);

    void multigridSubMatmulThread(const Range range, const MaterialGrid& materials, const std::vector<double>& vsub,
                            const std::vector<double>& vmul, std::vector<double>& vout);

    std::array<std::pair<int,float>,5> getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int i, int j);

    std::array<std::pair<int,float>,5> getMultigridMatrixEntriesForCell(const MaterialGrid &materials, int linearIdx);

    std::array<int,16> getFineGridStencilIdxs(const LinearIndexable2d &fineGridIndexer,int iCoarse, int jCoarse);

    std::array<int,4> getFineGridChildIdxs(const LinearIndexable2d &fineGridIndexer,int iCoarse, int jCoarse);

    std::array<float,4> getProlongationWeights(int finei, int finej);

    std::array<int,4> getCoarseProlongIdxs(const LinearIndexable2d &coarseGridIndexer, int finei, int finej);

    static const double m_tol;

    std::array<float,16> m_restrictionWeights;
    std::array<float,16> m_prolongationWeights;

    MaterialGrid& m_mainMaterialGrid;
    std::vector<MaterialGrid> m_materialSubgrids;
    std::vector<Grid2d<double>> m_pressureGrids;
    std::vector<Grid2d<double>> m_rhsGrids;

    void(*m_dampedJacobiThread)(LinearSolver*,const Range,const MaterialGrid&,std::vector<double>&,
                                 const std::vector<double>&, const std::vector<double>&);
};

#endif // PCGSOLVER_H
