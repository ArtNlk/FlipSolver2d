#include "grid2d.h"

#include <ios>
#include <sstream>

template<class T>
std::string Grid2d<T>::toString() const
{
    std::stringstream output;

    output << m_sizeI << 'x' << m_sizeJ << ':';

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            output << getAt(i,j) << ',';
        }
    }

    output.seekp(-1,std::ios::end);
    output << ';';

    return output.str();
}

template<>
std::string Grid2d<bool>::toString() const
{
    std::stringstream output;

    output << m_sizeI << 'x' << m_sizeJ << ':';

    for(int i = 0; i < m_sizeI; i++)
    {
        for(int j = 0; j < m_sizeJ; j++)
        {
            output << getAt(i,j) << ',';
        }
    }

    output.seekp(-1,std::ios::end);
    output << ';';

    return output.str();
}
