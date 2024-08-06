#ifndef INDEX2D_H
#define INDEX2D_H

#ifndef HAVE_SSIZE_T
#if defined(_WIN32) || defined(_WIN64)
typedef ptrdiff_t ssize_t; // ssize_t equivalent on Windows
#define HAVE_SSIZE_T
#elif defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/types.h>
#define HAVE_SSIZE_T
#else
// Fallback definition if no other platform-specific definition is available
typedef long ssize_t;
#define HAVE_SSIZE_T
#endif
#endif

struct Index2d
{
    Index2d() :
        i(0),
        j(0)
    {
    }

    Index2d(ssize_t i, ssize_t j) :
        i(i),
        j(j)
    {
    }

    ssize_t i;
    ssize_t j;
};

#endif // INDEX2D_H
