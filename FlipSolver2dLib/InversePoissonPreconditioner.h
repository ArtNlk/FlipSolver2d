#ifndef INVERSEPOISSONPRECONDITIONER_H
#define INVERSEPOISSONPRECONDITIONER_H

#include "linearindexable2d.h"
#include "matrixweights.h"
#include "threadpool.h"
#include "materialgrid.h"

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

    static InversePoissonPreconditioner create(double stepDt, double density, double dx, const MaterialGrid& materialGrid)
    {
        const double scale = stepDt / (density * dx * dx);

        InversePoissonPreconditioner output(materialGrid.linearSize()*0.33, materialGrid);

        const LinearIndexable2d& indexer = static_cast<LinearIndexable2d>(materialGrid);

        std::vector<Range> threadRanges = ThreadPool::i()->splitRange(materialGrid.linearSize());
        size_t currRangeIdx = 0;

        std::vector<std::array<double,2>> tempData(materialGrid.linearSize(),{0.0,0.0});

        for(size_t i = 0; i < materialGrid.sizeI(); i++)
        {
            for(size_t j = 0; j < materialGrid.sizeJ(); j++)
            {
                const size_t linIdx = indexer.linearIndex(i,j);

                if(!materialGrid.isFluid(i,j))
                {
                    tempData[linIdx] = {0.0,0.0};
                    continue;
                };

                double iDiag = materialGrid.inBounds(i-1, j) ? materialGrid.nonsolidNeighborCount(i-1, j) : 0.0;
                double jDiag = materialGrid.inBounds(i, j-1) ? materialGrid.nonsolidNeighborCount(i, j-1) : 0.0;
                iDiag = std::abs(iDiag) < 1e-9 ? 1 : iDiag * scale;
                jDiag = std::abs(jDiag) < 1e-9 ? 1 : jDiag * scale;
                double iNeg = materialGrid.inBounds(i-1, j) && materialGrid.isFluid(i-1, j) ? scale : 0.0;
                double jNeg = materialGrid.inBounds(i, j-1) && materialGrid.isFluid(i, j-1) ? scale : 0.0;

                tempData[linIdx] = {iNeg/iDiag, jNeg/jDiag};
            }
        }

        for(size_t i = 0; i < materialGrid.sizeI(); i++)
        {
            for(size_t j = 0; j < materialGrid.sizeJ(); j++)
            {
                const ssize_t linIdx = indexer.linearIndex(i,j);

                // if(!materialGrid.isFluid(i,j))
                // {
                //     if(linIdx >= threadRanges.at(currRangeIdx).end)
                //     {
                //         output.endThreadDataRange();
                //         currRangeIdx++;
                //     }
                //     continue;
                // }

                std::array<double,2> currRowData= tempData[linIdx];

                IndexedIPPCoefficientUnit unit;
                unit.unitIndex = linIdx;

                ssize_t b0Idx = linIdx - indexer.iLinearOffset();
                ssize_t b1Idx = linIdx - indexer.iLinearOffset() + indexer.jLinearOffset();
                ssize_t b2Idx = linIdx - indexer.jLinearOffset();
                ssize_t b3Idx = linIdx;
                ssize_t b4Idx = linIdx + indexer.jLinearOffset();
                ssize_t b5Idx = linIdx + indexer.iLinearOffset() - indexer.jLinearOffset();
                ssize_t b6Idx = linIdx + indexer.iLinearOffset();

                double b1data = tempData[b1Idx][1];
                double b4data = tempData[b4Idx][1];
                double b5data = tempData[b5Idx][0];
                double b6data = tempData[b6Idx][0];

                unit.data[0] = currRowData[0];
                unit.data[1] = currRowData[0] * b1data;
                unit.data[2] = currRowData[1];
                unit.data[3] = currRowData[0] * currRowData[0] + currRowData[1]*currRowData[1] + 1.0;
                unit.data[4] = b4data;
                unit.data[5] = currRowData[1] * b5data;
                unit.data[6] = b6data;

                unit.idx[0] = b0Idx;
                unit.idx[1] = b1Idx;
                unit.idx[2] = b2Idx;
                unit.idx[3] = b3Idx;
                unit.idx[4] = b4Idx;
                unit.idx[5] = b5Idx;
                unit.idx[6] = b6Idx;

                if(linIdx >= threadRanges.at(currRangeIdx).end)
                {
                    output.endThreadDataRange();
                    currRangeIdx++;
                }

                output.add(unit);
            }
        }

        if(currRangeIdx != threadRanges.size())
        {
            output.endThreadDataRange();
        }

        return output;
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
