#ifndef INVERSEPOISSONPRECONDITIONER_H
#define INVERSEPOISSONPRECONDITIONER_H

#include "linearindexable2d.h"
#include "matrixweights.h"
#include "threadpool.h"

#include <cstddef>
#include <limits>
#include <vector>

struct IndexedIPPCoefficientUnit
{
    IndexedIPPCoefficientUnit() :
        unitIndex(std::numeric_limits<size_t>::max()),
        data({0.0}),
        idx({-1})
    {

    }

    inline double multiply(const std::vector<double>& vec,
                           const LinearIndexable2d& indexer) const
    {
        double output = 0.0;
        for(int i = 0; i < data.size(); i++)
        {
            output += indexer.inBounds(idx[i]) ? data[i] * vec[idx[i]] : 0.0;
        }

        return output;
    }

    size_t unitIndex;
    std::array<double,7> data;
    std::array<ssize_t,7> idx;
};

class InversePoissonPreconditioner : public MatrixWeights<IndexedIPPCoefficientUnit>
{
public:
    InversePoissonPreconditioner(size_t size, const LinearIndexable2d& indexer) :
    MatrixWeights(indexer)
    {
        m_data.reserve(size);
    }

protected:
    void multiplyThread(Range vecRange, Range dataRange, const std::vector<double>& in, std::vector<double>& out) const
    {
        if(dataRange.size() == 0)
        {
            std::copy(in.begin() + vecRange.start, in.begin() + vecRange.end,out.begin() + vecRange.start);
            return;
        }

        //size_t nextVectorIndex = m_data[0].unitIndex;
        size_t nextDataIndex = dataRange.start;
        for(size_t idx = vecRange.start; idx < vecRange.end; idx++)
        {
            if(nextDataIndex < m_data.size() && nextDataIndex < dataRange.end
                && m_data[nextDataIndex].unitIndex == idx)
            {
                out[idx] = m_data[nextDataIndex].multiply(in, m_indexer);
                nextDataIndex++;
            }
            else
            {
                out[idx] = in[idx];
            }
        }
    }
};

#endif // INVERSEPOISSONPRECONDITIONER_H
