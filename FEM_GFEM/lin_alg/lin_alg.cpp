#include "lin_alg.h"

// resolver sistema linear K * U = F
// método de eliminação de Gauss
Vector Gauss_elimination(Matrix& K, Vector& F)
{
    std::size_t n {F.size()};
    Vector U(n);

    for (std::size_t i {0}; i < n; i++)
    {
        for (std::size_t j {0}; j < n; j++)
        {
            if (j != i)
            {
                double factor {K[j][i]/K[i][i]};
                for (std::size_t k {i}; k < n; k++)
                {
                    K[j][k] -= factor * K[i][k];
                }
                F[j] -= factor * F[i];
            }
        }
    }

    for (std::size_t i {0}; i < n; i++)
    {
        U[i] = F[i]/K[i][i];
    }
    return std::move(U);
}