#ifndef CUSTOMASSERT_H
#define CUSTOMASSERT_H

#include <iostream>
#include <vector>
#include <cmath>

#ifdef DEBUG
    #define ASSERT(cond)\
    {\
        if(!(cond)){\
            std::cout << "\nASSERT FAILED:"\
                        << #cond\
                        << " at: " << __FILE__ << " : "\
                        << __FUNCTION__ << ":"\
                        << __LINE__;\
            abort();\
        }\
    }
#else
    #define ASSERT(cond) ((void)0)
#endif

#ifdef DEBUG
    #define ASSERT_BETWEEN(val,min,max)\
        ASSERT((val) > (min) && (val) < (max))
#else
    #define ASSERT_BETWEEN(val,min,max) ((void)0)
#endif

#ifdef DEBUG
    template<typename T>
    bool isNanInf(T v)
    {
        static_assert(std::is_floating_point<T>::value);
        auto c = std::fpclassify(v);
        return (c == FP_NAN || c == FP_INFINITE);
    }
#else
    template<typename T>
    bool isNanInf(T v)
    {
        static_assert(std::is_floating_point<T>::value);
        return false;
    }
#endif

#ifdef DEBUG
    template<typename T>
    bool anyNanInf(std::vector<T> v)
    {
        static_assert(std::is_floating_point<T>::value);
        for(auto value : v)
        {
            if(isNanInf(value))
            {
                return true;
            }
        }
        return false;
    }
#else
    template<typename T>
    bool anyNanInf(std::vector<T> v)
    {
        static_assert(std::is_floating_point<T>::value);
        return false;
    }
#endif

#endif // CUSTOMASSERT_H
