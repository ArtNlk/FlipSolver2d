#ifndef LIGHTVISCOSITYWEIGHTS_H
#define LIGHTVISCOSITYWEIGHTS_H

#include "matrixweights.h"
#include "linearindexable2d.h"

struct LightViscosityWeightsUnit
{
    LightViscosityWeightsUnit() = default;

    inline double multiply(const std::vector<double>& vec, const LinearIndexable2d& indexer) const
    {
        const ssize_t im1Idx = indexer.linearIdxOfOffset(unitIndex,-1,0);
        const ssize_t ip1Idx = indexer.linearIdxOfOffset(unitIndex,1,0);
        const ssize_t jm1Idx = indexer.linearIdxOfOffset(unitIndex,0,-1);
        const ssize_t jp1Idx = indexer.linearIdxOfOffset(unitIndex,0,1);
        const ssize_t centerIdx = unitIndex;

        return diag*(centerIdx >= 0 && centerIdx < vec.size() ? vec.at(centerIdx) : 0.0) +
               im1*(im1Idx >= 0 && im1Idx < vec.size() ? vec.at(im1Idx) : 0.0) +
               ip1*(ip1Idx >= 0 && ip1Idx < vec.size() ? vec.at(ip1Idx) : 0.0) +
               jm1*(jm1Idx >= 0 && jm1Idx < vec.size() ? vec.at(jm1Idx) : 0.0) +
               jp1*(jp1Idx >= 0 && jp1Idx < vec.size() ? vec.at(jp1Idx) : 0.0);
    }

    size_t unitIndex;

    double diag;
    double im1;
    double ip1;
    double jm1;
    double jp1;
};

class LightViscosityWeights : public MatrixWeights<LightViscosityWeightsUnit>
{
public:
    LightViscosityWeights(LinearIndexable2d indexer):
    MatrixWeights(indexer)
    {}

protected:
    void multiplyThread(Range vecRange, Range dataRange, const std::vector<double>& in, std::vector<double>& out) const override
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

#endif // LIGHTVISCOSITYWEIGHTS_H
