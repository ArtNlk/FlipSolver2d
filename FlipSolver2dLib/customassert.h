#ifndef CUSTOMASSERT_H
#define CUSTOMASSERT_H

#include <iostream>

#ifdef DEBUG
    #define ASSERT(cond)\
    {\
        if(!cond){\
            std::cout << "ASERT FAILED:"\
                        << #cond\
                        << "at: " << __FILE__ << " : "\
                        << __FUNCTION__ << ":"\
                        << __LINE__;\
            abort();\
        }\
    }
#else
    #define ASSERT(cond) ((void)0)
#endif

#ifdef DEBUG
    #define ASSERT_BETWEEN_INC(val,min,max)\
        ASSERT(val > min && val < max)
#else
    #define ASSERT_BETWEEN(val,min,max) ((void)0)
#endif

#endif // CUSTOMASSERT_H
