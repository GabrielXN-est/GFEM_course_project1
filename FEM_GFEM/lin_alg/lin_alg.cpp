#include "lin_alg.h"
#include <thread>
#include <iostream>
#include <cmath>

void verify_if_singular(Matrix K)
{
    try
    {
        if (K.determinant() < std::pow(10, -12))
            throw std::runtime_error("Matrix is singular, cannot solve system");
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(1);
    }
}

// resolver sistema linear K * U = F
// método de eliminação de Gauss
Vector Gauss_elimination(Matrix K, Vector F)
{
    std::size_t n {F.size()};
    Vector U(n);

    std::thread t (verify_if_singular, K);

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

    t.join();

    for (std::size_t i {0}; i < n; i++)
    {
        U[i] = F[i]/K[i][i];
    }
    return std::move(U);
}