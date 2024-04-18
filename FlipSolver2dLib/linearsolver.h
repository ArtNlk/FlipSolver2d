#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <functional>
#include <unordered_map>
#include <vector>

#include "linearindexable2d.h"
#include "materialgrid.h"
#include "threadpool.h"
#include "staticmatrix.h"

class IPreconditioner
{
public:
    struct PreconditionerData
    {
        virtual ~PreconditionerData() = default;
    };

    virtual void apply(std::vector<double> const &in, std::vector<double> &out,
                       PreconditionerData* data = nullptr) = 0;
    virtual ~IPreconditioner() = default;
};

class StubPreconditioner : public IPreconditioner
{
public:
    void apply(std::vector<double> const &in, std::vector<double> &out,
               PreconditionerData* data = nullptr) override;
};

class IPPreconditioner : public IPreconditioner
{
public:
    IPPreconditioner(LinearIndexable2d indexer);

    void apply(std::vector<double> const &in, std::vector<double> &out,
               PreconditionerData* data = nullptr) override;

    struct IPPreconditionerData : public PreconditionerData
    {
        IPPreconditionerData(StaticMatrix &m) :
            m(m)
        {

        }

        StaticMatrix &m;
    };

private:
    void firstStepIPPMatmulThread(Range r, const std::vector<double> &in,
                                  std::vector<double> &out, StaticMatrix& m);

    void secondStepIPPMatmulThread(Range r, const std::vector<double> &in,
                                   std::vector<double> &out, StaticMatrix& m);

    LinearIndexable2d m_indexer;
};

class LinearSolver
{
public:
    LinearSolver(MaterialGrid& materialGrid, int maxMultigridDepth);
    using SparseMatRowElements = std::array<std::pair<int,double>,5>;

    bool solve(const StaticMatrix &matrixIn, std::vector<double> &result,
               const std::vector<double> &vec, std::shared_ptr<IPreconditioner> precond, IPreconditioner::PreconditionerData* data,
               int iterLimit = 20, double tol = 1e-6);

    friend class LinearSolver_sse42;

protected:
    void applyICPrecond(const DynamicMatrix &precond, std::vector<double> const &in, std::vector<double> &out);
    DynamicMatrix calcPrecond(const DynamicMatrix &matrix);
};

#endif // PCGSOLVER_H
