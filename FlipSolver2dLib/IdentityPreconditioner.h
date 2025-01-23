#ifndef IDENTITYPRECONDITIONER_H
#define IDENTITYPRECONDITIONER_H

#include "matrixweights.h"

struct IdentityPreconditionerUnit
{

};

class IdentityPreconditioner : public MatrixWeights<IdentityPreconditionerUnit>
{
public:
    IdentityPreconditioner() :
    MatrixWeights(LinearIndexable2d(0,0))
    {

    }

    void multiply(const std::vector<double>& in, std::vector<double>& out) const override
    {
        out = in;
    }

protected:
    void multiplyThread(Range vecRange, Range dataRange, const std::vector<double>& in, std::vector<double>& out) const override
    {
        return;
    }
};

#endif // IDENTITYPRECONDITIONER_H
