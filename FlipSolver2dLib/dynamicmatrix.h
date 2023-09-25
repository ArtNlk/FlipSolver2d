#ifndef DYNAMICMATRIX_H
#define DYNAMICMATRIX_H

#include <vector>
#include <utility>
#include <sstream>

#include "linearindexable2d.h"

class Logger;

using SparseRowUnit = std::pair<int, double>;
using SparseRow = std::vector<SparseRowUnit>;

class DynamicMatrix
{
public:
    friend Logger;
    
    DynamicMatrix(int size, int avgRowLength = 7);
    
    DynamicMatrix(const DynamicMatrix &m) = default;
    
    DynamicMatrix(DynamicMatrix &&m) = default;
    
    void copyUpperTriangleTo(DynamicMatrix &m) const;

    void setAdiag(int i, int j, double value, LinearIndexable2d &indexer);

    void setAx(int i, int j, double value, LinearIndexable2d &indexer);

    void setAy(int i, int j, double value, LinearIndexable2d &indexer);

    void addTo(int i, int j, double value);

    void addToAdiag(int i, int j, double value, LinearIndexable2d &indexer);

    void addToAx(int i, int j, double value, LinearIndexable2d &indexer);

    void addToAy(int i, int j, double value, LinearIndexable2d &indexer);

    int rowSize(int rowIndex) const;

    bool isStored(int rowIdx, int colIdx);

    int elementCount() const ;

    std::vector<SparseRow>& data();

    void setValue(int rowIndex, int columnIndex, double value);
    double getValue(int rowIndex, int columnIndex) const;

    std::vector<double> operator*(std::vector<double> &v) const;

    std::string toString();

    size_t size() const;

protected:
    std::vector<SparseRow> m_rows;
    size_t m_size;
    int m_elementCount;
};

#endif // DYNAMICMATRIX_H
