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
        iNeg(0.0),
        iPos(0.0),
        jNeg(0.0),
        jPos(0.0)
    {

    }

    inline double multiply(const std::vector<double>& vec,
                           const LinearIndexable2d& indexer) const
    {
        const ssize_t jNegLinIdx = indexer.linearIdxOfOffset(unitIndex,0,-1);
        const ssize_t jPosLinIdx = indexer.linearIdxOfOffset(unitIndex,0,1);
        const ssize_t iNegLinIdx = indexer.linearIdxOfOffset(unitIndex,-1,0);
        const ssize_t iPosLinIdx = indexer.linearIdxOfOffset(unitIndex,1,0);
        const ssize_t centerIdx = unitIndex;

        const double diag = 1.0 + iNeg*iNeg + jNeg*jNeg;

        return  diag * (centerIdx >= 0 && centerIdx < vec.size() ? vec.at(centerIdx) : 0.0) +
                jNeg * (jNegLinIdx >= 0 && jNegLinIdx < vec.size() ? vec.at(jNegLinIdx) : 0.0) +
                jPos * (jPosLinIdx >= 0 && jPosLinIdx < vec.size() ? vec.at(jPosLinIdx) : 0.0) +
                iNeg * (iNegLinIdx >= 0 && iNegLinIdx < vec.size() ? vec.at(iNegLinIdx) : 0.0) +
                iPos * (iPosLinIdx >= 0 && iPosLinIdx < vec.size() ? vec.at(iPosLinIdx) : 0.0);

    }

    size_t unitIndex;
    double iNeg;
    double iPos;
    double jNeg;
    double jPos;
};

class InversePoissonPreconditioner : public MatrixWeights<IndexedIPPCoefficientUnit>
{
public:
    InversePoissonPreconditioner(size_t size, LinearIndexable2d& indexer) :
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
