#include "lin_alg.h"
 // estudar sintaxe
Matrix::Matrix (Vector v) : mat {((!v.transposed)? 
    std::vector<std::vector<double>>(v.size(), std::vector<double>(1, 0.))
    : std::vector<std::vector<double>>(1, std::vector<double>(v.size(), 0.)))}
{
    if (!v.transposed)
        for (std::size_t i {0}; i < v.size(); i++)
        {
            mat[i][0] = v.get(i);
        }
    else
        for (std::size_t i {0}; i < v.size(); i++)
        {
            mat[0][i] = v.get(i);
        }
}

Matrix operator* (double n, Matrix&& m)
{
    for (std::vector<double>& row : m.mat)
    {
        for (double& val : row)
        {
            val *= n;
        }
    }
    return std::move(m);
}