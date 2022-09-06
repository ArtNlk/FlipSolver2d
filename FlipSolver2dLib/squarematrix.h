#ifndef IMATRIX_H
#define IMATRIX_H

#include "linearindexable2d.h"

class SquareMatrix : public LinearIndexable2d
{
public:
    SquareMatrix(int size) :
    LinearIndexable2d(size,size),
    m_elementCount(0)
    {}

    virtual double getValue(int row, int column) const = 0;
    virtual void setValue(int row, int column, double value) = 0;
    int size() const {return m_sizeI;};
    virtual int rowSize(int rowIndex) = 0;
    int elementCount() {return m_elementCount;};

protected:
    int m_elementCount;
};

#endif // IMATRIX_H
