#ifndef UPPERTRIANGULARMATRIX_H
#define UPPERTRIANGULARMATRIX_H

#include "sparsematrix.h"
#include "dynamicuppertriangularsparsematrix.h"

class UpperTriangularMatrix : public SparseMatrix
{
public:
    UpperTriangularMatrix(const DynamicUpperTriangularSparseMatrix &dynamicMatrix);

    virtual double getValue(int row, int col) const override;
};

#endif // UPPERTRIANGULARMATRIX_H
