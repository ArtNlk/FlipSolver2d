#include "dynamicmatrix.h"

#include <algorithm>
#include <optional>

#include "linearindexable2d.h"
#include "mathfuncs.h"

DynamicMatrix::DynamicMatrix(int size, int avgRowLength) :
    m_rows(size),
    m_size(size),
    m_elementCount(0)
{
    for(int i = 0; i < size; i++)
    {
        m_rows[i].reserve(avgRowLength);
    }
}

void DynamicMatrix::copyUpperTriangleTo(DynamicMatrix &m) const
{
    ASSERT(m.size() == this->size());
    int size = m.size();
    m.m_elementCount = 0;

    for(int rowIdx = 0; rowIdx < size; rowIdx++)
    {
        for(int rowElementIdx = 0; rowElementIdx < m_rows[rowIdx].size(); rowElementIdx++)
        {
            if(m_rows[rowIdx][rowElementIdx].first >= rowIdx)
            {
                int leftoverLength = m_rows[rowIdx].size() - rowElementIdx;
                auto last = m_rows[rowIdx].cend();
                auto first = last - leftoverLength;
                m.m_rows[rowIdx] = std::vector(first,last);
                m.m_elementCount += m.m_rows[rowIdx].size();
            }
        }
    }
}


void DynamicMatrix::setAdiag(int i, int j, double value, LinearIndexable2d &indexer)
{
    ASSERT_BETWEEN(i,-2,m_size);
    ASSERT_BETWEEN(j,-2,m_size);
    int index = indexer.linearIndex(i,j);
    setValue(index, index, value);
}

void DynamicMatrix::setAx(int i, int j, double value, LinearIndexable2d &indexer)
{
    int rowIndex = indexer.linearIndex(i,j);
    int colIndex = indexer.linearIndex(i+1,j);
    if(colIndex == -1 || rowIndex == -1)
    {
        return;
    }

    setValue(rowIndex,colIndex, value);
    setValue(colIndex,rowIndex, value);
}

void DynamicMatrix::setAy(int i, int j, double value, LinearIndexable2d &indexer)
{
    int rowIndex = indexer.linearIndex(i,j);
    int colIndex = indexer.linearIndex(i,j+1);
    if(colIndex == -1 || rowIndex == -1)
    {
        return;
    }

    setValue(rowIndex,colIndex, value);
    setValue(colIndex,rowIndex, value);
}

void DynamicMatrix::addTo(int i, int j, double value)
{
    ASSERT_BETWEEN(i,-2,m_size);
    ASSERT_BETWEEN(j,-2,m_size);
    double temp = getValue(i,j);
    setValue(i,j, temp + value);
   // setValue(j,i, temp + value);
}

void DynamicMatrix::addToAdiag(int i, int j, double value, LinearIndexable2d &indexer)
{
    ASSERT_BETWEEN(i,-2,m_size);
    ASSERT_BETWEEN(j,-2,m_size);
    int index = indexer.linearIndex(i,j);
    setValue(index, index, getValue(index,index) + value);
}

void DynamicMatrix::addToAx(int i, int j, double value, LinearIndexable2d &indexer)
{
    ASSERT_BETWEEN(i,-2,m_size);
    ASSERT_BETWEEN(j,-2,m_size);
    int rowIndex = indexer.linearIndex(i,j);
    int colIndex = indexer.linearIndex(i+1,j);

    setValue(rowIndex,colIndex, getValue(rowIndex, colIndex) + value);
}

void DynamicMatrix::addToAy(int i, int j, double value, LinearIndexable2d &indexer)
{
    ASSERT_BETWEEN(i,-2,m_size);
    ASSERT_BETWEEN(j,-2,m_size);
    int rowIndex = indexer.linearIndex(i,j);
    int colIndex = indexer.linearIndex(i,j+1);

    return setValue(rowIndex,colIndex, getValue(rowIndex, colIndex) + value);
}

int DynamicMatrix::rowSize(int rowIndex) const { return m_rows[rowIndex].size();}

bool DynamicMatrix::isStored(int rowIdx, int colIdx)
{
    ASSERT_BETWEEN(rowIdx, -1, m_size);
    ASSERT_BETWEEN(colIdx, -1, m_size);

    auto sparseRowUnitCmp = [](const SparseRowUnit &a, const SparseRowUnit &b){
        return a.first < b.first;
    };

    std::vector<SparseRowUnit> selectedRow = m_rows[rowIdx];

    SparseRowUnit searchValue(colIdx,0.0);

    return std::binary_search(selectedRow.begin(),selectedRow.end(),searchValue,sparseRowUnitCmp);
}

int DynamicMatrix::elementCount() const { return m_elementCount;}

std::vector<SparseRow>& DynamicMatrix::data() { return m_rows;}

void DynamicMatrix::setValue(int rowIndex, int columnIndex, double value)
{
    SparseRow &targetRow = m_rows[rowIndex];
    for(int i = 0; i < targetRow.size(); i++)
    {
        if(targetRow[i].first == columnIndex)
        {
            targetRow[i].second = value;
            return;
        }
        if(targetRow[i].first > columnIndex)
        {
            targetRow.insert(targetRow.begin()+i,SparseRowUnit(columnIndex,value));
            m_elementCount++;
            return;
        }
    }
    targetRow.push_back(SparseRowUnit(columnIndex,value));
    m_elementCount++;
}

double DynamicMatrix::getValue(int rowIndex, int columnIndex) const
{
    const SparseRow &targetRow = m_rows[rowIndex];
    for(const auto & column : targetRow)
    {
        if(column.first == columnIndex)
        {
            return column.second;
        }
    }

    return 0;
}

std::vector<double> DynamicMatrix::operator*(std::vector<double> &v) const
{
    std::vector<double> output(v.size());
    int rowIdx = 0;
    for(const auto &row : m_rows)
    {
        double result = 0.0;
        for(const auto &pair : row)
        {
            result += v[pair.first] * pair.second;
        }
        output[rowIdx] = result;
        rowIdx++;
    }

    return output;
}

std::string DynamicMatrix::toString()
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

size_t DynamicMatrix::size() const
{
    return m_size;
}
