#ifndef NPOINTSTENCILMATRIX_H
#define NPOINTSTENCILMATRIX_H

#include <type_traits>
#include <vector>

template<size_t NPoints>
class NPointStencilMatrix
{
public:
    NPointStencilMatrix(size_t size) :
        m_vals(size * NPoints, 0.0),
        m_idxs(size * NPoints, 0.0)
    {
    };

    template<size_t np = NPoints, std::enable_if_t<NPoints == 5>>
    void vmul(const std::vector<double>& vin, std::vector<double>& vout)
    {
        for(size_t idx = 0; idx < vin.size(); idx++)
        {
            double temp = 0.0;
            temp += vin[m_idxs[idx*m_stencilSize + 0]] * m_vals[idx*m_stencilSize + 0];
            temp += vin[m_idxs[idx*m_stencilSize + 1]] * m_vals[idx*m_stencilSize + 1];
            temp += vin[m_idxs[idx*m_stencilSize + 2]] * m_vals[idx*m_stencilSize + 2];
            temp += vin[m_idxs[idx*m_stencilSize + 3]] * m_vals[idx*m_stencilSize + 3];
            temp += vin[m_idxs[idx*m_stencilSize + 4]] * m_vals[idx*m_stencilSize + 4];
            vout[idx] = temp;
        }
    }

    template<size_t np = NPoints, std::enable_if_t<NPoints != 5>>
    void vmul(const std::vector<double>& vin, std::vector<double>& vout)
    {
        for(size_t idx = 0; idx < vin.size(); idx++)
        {
            double temp = 0.0;
            for(size_t pIdx = 0; pIdx < m_stencilSize; pIdx++)
            {
                temp += vin[m_idxs[idx*m_stencilSize + pIdx]] * m_vals[idx*m_stencilSize + pIdx];
            }
            vout[idx] = temp;
        }
    }

    void setStencil(size_t stencilIdx)
    {

    }

    std::vector<double>& values();
    std::vector<size_t>& indexes();

protected:
    std::vector<double> m_vals;
    std::vector<size_t> m_idxs;
    const size_t m_stencilSize = NPoints;
};

#endif // NPOINTSTENCILMATRIX_H
