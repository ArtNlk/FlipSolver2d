#ifndef STATICMATRIX_H
#define STATICMATRIX_H

#include "dynamicmatrix.h"
#include "matrixbase.h"

#include "threadpool.h"

class StaticMatrix : public MatrixBase
{
public:
    template<size_t MaxRowSize>
    static StaticMatrix fromDynamic(const DynamicMatrix<MaxRowSize> &in)
    {
        using DataUnit = DynamicMatrix<MaxRowSize>::SparseRowDataUnit;
        using IndexUnit = DynamicMatrix<MaxRowSize>::SparseRowIndexUnit;
        const size_t reserveSize = MaxRowSize * in.size() / 2;

        StaticMatrix output;
        output.m_size = in.size();
        output.m_indexes.reserve(reserveSize);
        output.m_values.reserve(reserveSize);
        output.m_rowStart.reserve(in.size() + 1);

        for (size_t rowIdx = 0; rowIdx < in.size(); rowIdx++) {
            const DataUnit& rowDataUnit = in.data()[rowIdx];
            const IndexUnit& rowIndexUnit = in.indexes()[rowIdx];

            output.m_rowStart.push_back(output.m_values.size());

            if (rowDataUnit.isEmpty())
            {
                continue;
            }

            for(size_t rowElementIdx = 0; rowElementIdx < rowDataUnit.size(); rowElementIdx++)
            {
                output.m_indexes.push_back(rowIndexUnit.data()[rowElementIdx]);
                output.m_values.push_back(rowDataUnit.data()[rowElementIdx]);
            }
        }

        output.m_rowStart.push_back(output.m_values.size());

        return output;
    }

    void multiply(const std::vector<double>& in, std::vector<double>& out) const override;

    std::string toString();

    size_t size() const;

protected:
    StaticMatrix();

    void mulThread(Range range, const std::vector<double> &vin, std::vector<double>& vout) const;

    std::vector<size_t> m_indexes;
    std::vector<double> m_values;
    std::vector<size_t> m_rowStart;

    size_t m_size;
};

#endif // STATICMATRIX_H
