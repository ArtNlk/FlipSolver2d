#ifndef INDEX2D_H
#define INDEX2D_H


struct Index2d
{
    Index2d() :
        i(0),
        j(0)
    {
    }

    Index2d(int i, int j) :
        i(i),
        j(j)
    {
    }

    int i;
    int j;
};

#endif // INDEX2D_H
