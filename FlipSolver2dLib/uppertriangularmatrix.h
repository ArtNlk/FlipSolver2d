#ifndef UPPERTRIANGULARMATRIX_H
#define UPPERTRIANGULARMATRIX_H

#include "dynamicuppertriangularsparsematrix.h"
#include "linearindexable2d.h"

#include "threadpool.h"

using StaticRowUnit = std::pair<int, double>;

class UpperTriangularMatrix : public SquareMatrix
{
public:

    UpperTriangularMatrix(const DynamicUpperTriangularSparseMatrix &dynamicMatrix);

    double getValue(int row, int col) const override;

    void setValue(int row, int col, double value) override;

    int rowSize(int rowIndex) const override;

    double Adiag(int i, int j, LinearIndexable2d &indexer) const;

    double Ax(int i, int j, LinearIndexable2d &indexer) const;

    double Ay(int i, int j, LinearIndexable2d &indexer) const;

    std::vector<double> operator*(std::vector<double> &v) const;

    std::string toString();

    void mulThread(Range range, std::vector<double> &vin, std::vector<double>& vout) const;
protected:
    std::vector<StaticRowUnit> m_values;
    std::vector<int> m_rowStart;
};

#endif // UPPERTRIANGULARMATRIX_H
