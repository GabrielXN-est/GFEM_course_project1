#include "lin_alg.h"

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