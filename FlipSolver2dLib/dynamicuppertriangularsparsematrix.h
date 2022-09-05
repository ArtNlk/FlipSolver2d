#ifndef DYNAMICUPPERTRIANGULARSPARSEMATRIX_H
#define DYNAMICUPPERTRIANGULARSPARSEMATRIX_H

#include <vector>
#include <utility>
#include <sstream>

#include "fluidgrid.h"
#include "squarematrix.h"

class Logger;

class DynamicUpperTriangularSparseMatrix : public SquareMatrix
{
public:
    typedef std::pair<int,double> SparseRowUnit;
    typedef std::vector<SparseRowUnit> SparseRow;
    friend Logger;

    DynamicUpperTriangularSparseMatrix(int size, int avgRowLength = 7);

    void setAdiag(int i, int j, double value, MACFluidGrid &grid);

    void setAx(int i, int j, double value, MACFluidGrid &grid);

    void setAy(int i, int j, double value, MACFluidGrid &grid);

    void addTo(int i, int j, double value);

    void addToAdiag(int i, int j, double value, MACFluidGrid &grid);

    void addToAx(int i, int j, double value, MACFluidGrid &grid);

    void addToAy(int i, int j, double value, MACFluidGrid &grid);

    virtual int rowSize(int rowIndex) override;

    int elementCount() ;

    const std::vector<SparseRow> data() const;

    virtual void setValue(int rowIndex, int columnIndex, double value) override;
    virtual double getValue(int rowIndex, int columnIndex) const override;

    std::string toString();

protected:

    void internalSetValue(int rowIndex, int columnIndex, double value);

    std::vector<SparseRow> m_rows;
    int m_size;
    int m_elementCount;
};

#endif // DYNAMICUPPERTRIANGULARSPARSEMATRIX_H
