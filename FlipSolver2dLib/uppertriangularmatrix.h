#ifndef UPPERTRIANGULARMATRIX_H
#define UPPERTRIANGULARMATRIX_H

#include "dynamicuppertriangularsparsematrix.h"

class UpperTriangularMatrix : public SquareMatrix
{
public:
    typedef std::pair<int,double> StaticRowUnit;

    UpperTriangularMatrix(DynamicUpperTriangularSparseMatrix &dynamicMatrix);

    virtual double getValue(int row, int col) const override;

    virtual void setValue(int row, int col, double value) override;

    virtual int rowSize(int rowIndex) override;

    double Adiag(int i, int j, MACFluidGrid &grid) const;

    double Ax(int i, int j, MACFluidGrid &grid) const;

    double Ay(int i, int j, MACFluidGrid &grid) const;

    std::vector<double> operator*(  std::vector<double> &v) const;

    std::string toString();

protected:
    std::vector<StaticRowUnit> m_values;
    std::vector<int> m_rowStart;
    int m_size;
};

#endif // UPPERTRIANGULARMATRIX_H
