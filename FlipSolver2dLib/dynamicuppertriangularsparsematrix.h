#ifndef DYNAMICUPPERTRIANGULARSPARSEMATRIX_H
#define DYNAMICUPPERTRIANGULARSPARSEMATRIX_H

#include <vector>
#include <utility>
#include <sstream>

#include "squarematrix.h"

class Logger;

using SparseRowUnit = std::pair<int, double>;
using SparseRow = std::vector<SparseRowUnit>;

class DynamicUpperTriangularSparseMatrix : public SquareMatrix
{
public:
    friend Logger;

    DynamicUpperTriangularSparseMatrix(int size, int avgRowLength = 7);

    DynamicUpperTriangularSparseMatrix(const DynamicUpperTriangularSparseMatrix &m) = default;

    DynamicUpperTriangularSparseMatrix(DynamicUpperTriangularSparseMatrix &&m) = default;

    void copyUpperTriangleTo(DynamicUpperTriangularSparseMatrix &m) const;

    void setAdiag(int i, int j, double value, LinearIndexable2d &indexer);

    void setAx(int i, int j, double value, LinearIndexable2d &indexer);

    void setAy(int i, int j, double value, LinearIndexable2d &indexer);

    void addTo(int i, int j, double value);

    void addToAdiag(int i, int j, double value, LinearIndexable2d &indexer);

    void addToAx(int i, int j, double value, LinearIndexable2d &indexer);

    void addToAy(int i, int j, double value, LinearIndexable2d &indexer);

    int rowSize(int rowIndex) const override;

    bool isStored(int rowIdx, int colIdx);

    int elementCount() const ;

    std::vector<SparseRow>& data();

    void setValue(int rowIndex, int columnIndex, double value) override;
    double getValue(int rowIndex, int columnIndex) const override;

    std::vector<double> operator*(std::vector<double> &v) const;

    std::string toString();

protected:
    std::vector<SparseRow> m_rows;
    int m_size;
    int m_elementCount;
};

#endif // DYNAMICUPPERTRIANGULARSPARSEMATRIX_H
