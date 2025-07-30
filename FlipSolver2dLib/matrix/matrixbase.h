#ifndef MATRIXBASE_H
#define MATRIXBASE_H

#include <vector>

class MatrixBase
{
public:
    MatrixBase() = default;
    virtual ~MatrixBase() = default;
    virtual void multiply(const std::vector<double>& in, std::vector<double>& out) const = 0;
};

#endif // MATRIXBASE_H
