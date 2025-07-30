#ifndef DYNAMICMATRIX_H
#define DYNAMICMATRIX_H

#include <vector>
#include <string>
#include <sstream>

#include "customassert.h"

class Logger;

template<size_t MaxRowSize, class StorageType>
class SparseRowUnit
{
public:
    SparseRowUnit():
        m_currentElementCount(0)
    {

    }

    size_t addElement(StorageType value)
    {
        ASSERT(m_currentElementCount+1 < m_array.size());

        size_t newIdx = m_currentElementCount;

        m_data[m_currentElementCount] = value;
        m_currentElementCount++;
        return newIdx;
    }

    bool isEmpty() const
    {
        return m_currentElementCount == 0;
    }

    size_t size() const
    {
        return m_currentElementCount;
    }

    size_t maxSize() const
    {
        return m_data.size();
    }

    StorageType& at(size_t arrayIndex)
    {
        ASSERT(arrayIndex < m_data.size());

        return m_data[arrayIndex];
    }

    const std::array<SparseRowUnit, MaxRowSize>& data() const
    {
        return m_data;
    }

private:
    size_t m_currentElementCount;
    std::array<SparseRowUnit, MaxRowSize> m_data;
};

template<size_t MaxRowSize>
class DynamicMatrix
{
public:
    friend Logger;

    using SparseRowIndexUnit = SparseRowUnit<MaxRowSize, size_t>;
    using SparseRowDataUnit = SparseRowUnit<MaxRowSize, double>;
    
    DynamicMatrix(size_t size):
        m_indexes(size),
        m_data(size),
        m_elementCount(0)
    {

    }

    size_t rowSize(size_t rowIndex) const
    {
        return m_data[rowIndex].size();
    }

    void addValue(size_t rowIndex, size_t columnIndex, double value)
    {
        m_indexes[rowIndex].addElement(columnIndex);
        m_data[rowIndex].addElement(value);
    }

    void addTo(size_t rowIndex, size_t columnIndex, double delta)
    {
        for(size_t storeIndex = 0; storeIndex < m_indexes[rowIndex]; storeIndex++)
        {
            if(m_indexes[rowIndex].data()[storeIndex] == columnIndex)
            {
                m_data[rowIndex].data()[storeIndex] += delta;
                break;
            }
        }
    }

    std::string toString()
    {
        std::ostringstream output;
        const size_t size = m_data.size();
        for(size_t row = 0; row < size; row++)
        {
            for(size_t elementIndex = 0; elementIndex < m_data[row].size(); elementIndex++)
            {
                output << row << ' ' << m_indexes[row][elementIndex] << ' ' << m_data[row][elementIndex] << '\n';
            }
        }

        return output.str();
    }

    size_t size() const
    {
        return m_data.size();
    }

    const std::vector<SparseRowIndexUnit>& indexes() const
    {
        return m_indexes;
    }

    const std::vector<SparseRowDataUnit>& data() const
    {
        return m_data;
    }

protected:
    std::vector<SparseRowIndexUnit> m_indexes;
    std::vector<SparseRowDataUnit> m_data;
    int m_elementCount;
};

#endif // DYNAMICMATRIX_H
