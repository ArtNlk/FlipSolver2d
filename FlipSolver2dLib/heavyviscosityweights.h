#ifndef HEAVYVISCOSITYWEIGHTS_H
#define HEAVYVISCOSITYWEIGHTS_H

#include "linearindexable2d.h"
#include "matrixweights.h"
#include <array>

struct HeavyViscosityWeightsUnit
{
    HeavyViscosityWeightsUnit() = default;

    inline double multiply(const std::vector<double>& vec, const LinearIndexable2d& indexer) const
    {
        return 0;
    }

    size_t unitIndex;
    std::array<double, 9> data;
    std::array<ssize_t, 9> indexes;
};

class HeavyViscosityWeights : public MatrixWeights<HeavyViscosityWeightsUnit>
{
public:
    HeavyViscosityWeights(size_t reserveSize) :
        MatrixWeights<HeavyViscosityWeightsUnit>(LinearIndexable2d(0,0)) //Dummy
    {
        m_data.reserve(reserveSize);
    };

    void multiply(const std::vector<double> &in, std::vector<double> &out) const override;

    void multiplyThread(Range vecRange, Range dataRange, const std::vector<double>& in, std::vector<double>& out) const override
    {
        for(size_t idx = vecRange.start; idx < vecRange.end; idx++)
        {
            if(m_data[idx].unitIndex == idx)
            {
                out[idx] = m_data[idx].multiply(in, m_indexer);
            }
            else
            {
                out[idx] = in[idx];
            }
        }
    }
};

#endif // HEAVYVISCOSITYWEIGHTS_H
