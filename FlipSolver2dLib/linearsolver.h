#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <functional>
#include <unordered_map>
#include <vector>

#include "linearindexable2d.h"
#include "materialgrid.h"
#include "threadpool.h"
#include "uppertriangularmatrix.h"

class LinearSolver
{
public:
    LinearSolver(MaterialGrid& materialGrid, int maxMultigridDepth);
    using SparseMatRowElements = std::array<std::pair<int,double>,5>;
    using MatElementProvider = std::function<SparseMatRowElements(int)>;

    bool solve(const DynamicUpperTriangularSparseMatrix &matrixIn, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 20);
    bool mfcgSolve(MatElementProvider elementProvider, std::vector<double> &result, const std::vector<double> &vec, int iterLimit = 20);

protected:
    void nomatVMul(MatElementProvider elementProvider, const std::vector<double> &vin, std::vector<double> &vout);
    void nomatVMulThread(Range range, MatElementProvider elementProvider,
                         const std::vector<double> &vin, std::vector<double> &vout);
    void applyMGPrecond(std::vector<double> const &in, std::vector<double> &out);
    void applyICPrecond(const DynamicUpperTriangularSparseMatrix &precond, std::vector<double> const &in, std::vector<double> &out);
    DynamicUpperTriangularSparseMatrix calcPrecond(const DynamicUpperTriangularSparseMatrix &matrix);

    void updateSubgrids();

    void propagateMaterialGrid(const MaterialGrid& fineGrid, MaterialGrid &coarseGrid);

    void restrictGrid(const MaterialGrid& coarseMaterials, const Grid2d<float> &fineGrid, Grid2d<float> coarseGrid);

    void prolongateGrid(const MaterialGrid& fineMaterials, const Grid2d<float> &coarseGrid, Grid2d<float> fineGrid);

    std::array<int,16> getFineGridStencilIdxs(const LinearIndexable2d &fineGridIndexer,int iCoarse, int jCoarse);

    std::array<int,4> getFineGridChildIdxs(const LinearIndexable2d &fineGridIndexer,int iCoarse, int jCoarse);

    std::array<float,4> getProlongationWeights(int finei, int finej);

    std::array<int,4> getCoarseProlongIdxs(const LinearIndexable2d &coarseGridIndexer, int finei, int finej);

    static const double m_tol;

    std::array<float,16> m_restrictionWeights;
    std::array<float,16> m_prolongationWeights;

    MaterialGrid& m_mainMaterialGrid;
    std::vector<MaterialGrid> m_materialSubgrids;
    std::vector<Grid2d<float>> m_pressureGrids;
};

#endif // PCGSOLVER_H
