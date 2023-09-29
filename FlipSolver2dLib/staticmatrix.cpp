#include "staticmatrix.h"
#include "logger.h"

#include "threadpool.h"

#include <stdexcept>
#include <vector>

StaticMatrix::StaticMatrix(const DynamicMatrix &dynamicMatrix):
    m_rowStart(dynamicMatrix.size()+1),
    m_size(dynamicMatrix.size())
{
    m_values.resize(dynamicMatrix.elementCount());

    m_rowStart[0] = 0;
    m_rowStart[size()] = dynamicMatrix.elementCount();
    for(int i = 1; i < size(); i++)
    {
        m_rowStart[i] = m_rowStart[i-1] + dynamicMatrix.rowSize(i-1);
    }

    int totalIndex = 0;
    std::vector<SparseRow> rows = const_cast<DynamicMatrix&>(dynamicMatrix).data();
    for(int i = 0; i < size(); i++)
    {
        for(int j = 0; j < rows[i].size(); j++)
        {
            m_values[totalIndex] = rows[i][j];
            totalIndex++;
        }
    }
}

double StaticMatrix::getValue(int row, int col) const
{
    if(row == -1 || col == -1 || !(row<m_size && row>=0 && col<m_size && col >=0)) {return 0.0;}
    int rowStartIndex = m_rowStart[row];
    int rowEndIndex = m_rowStart[row+1];
    for(int i = rowStartIndex; i < rowEndIndex; i++)
    {
        if(m_values[i].first == col)
        {
            return m_values[i].second;
        }
    }

    return 0;
}

double StaticMatrix::Adiag(int i, int j, LinearIndexable2d &indexer) const
{
    ASSERT_BETWEEN(i,-2,m_size);
    ASSERT_BETWEEN(j,-2,m_size);
    if(i < 0 || j < 0)
    {
        return 0.0;
    }
    int index = indexer.linearIndex(i,j);
    return getValue(index,index);
}

double StaticMatrix::Ax(int i, int j, LinearIndexable2d &indexer) const
{
    ASSERT_BETWEEN(i,-2,m_size);
    ASSERT_BETWEEN(j,-2,m_size);
    if(i < 0 || j < 0)
    {
        return 0.0;
    }
    int rowIndex = indexer.linearIndex(i,j);
    int colIndex = indexer.linearIndex(i+1,j);

    return getValue(rowIndex,colIndex);
}

double StaticMatrix::Ay(int i, int j, LinearIndexable2d &indexer) const
{
    ASSERT_BETWEEN(i,-2,m_size);
    ASSERT_BETWEEN(j,-2,m_size);
    if(i < 0 || j < 0)
    {
        return 0.0;
    }
    int rowIndex = indexer.linearIndex(i,j);
    int colIndex = indexer.linearIndex(i,j+1);

    return getValue(rowIndex,colIndex);
}

std::vector<double> StaticMatrix::operator*(std::vector<double> &v) const
{
    std::vector<double> output(v.size(),0.0);

//    for(int i = 0; i < size(); i++)
//    {
//        //vout[i] = 0;
//        for(int j = m_rowStart[i]; j < m_rowStart[i+1]; j++)
//        {
//            output[i] += m_values[j].second * v[m_values[j].first];
//        }
//    }

    std::vector<Range> ranges = ThreadPool::i()->splitRange(output.size());

    for(const Range& range : ranges)
    {
        ThreadPool::i()->enqueue(&StaticMatrix::mulThread,this,range,
                                 std::ref(v),std::ref(output));
    }
    ThreadPool::i()->wait();

    return output;
}

void StaticMatrix::mulThread(Range range, std::vector<double>& vin, std::vector<double> &vout) const
{
    for(int i = range.start; i < range.end; i++)
    {
        //vout[i] = 0;
        for(int j = m_rowStart[i]; j < m_rowStart[i+1]; j++)
        {
            vout[i] += m_values[j].second * vin[m_values[j].first];
        }
    }
}

std::string StaticMatrix::toString()
{
    std::ostringstream output;
    for(int i = 0; i < m_size*m_size; i++)
    {
        output << "|";
        for(int j = 0; j < m_size*m_size; j++)
        {
            output << "\t" << getValue(i,j) << ",";
        }
        output << "|\n";
    }

    return output.str();
}

size_t StaticMatrix::size() const
{
    return m_size;
}
