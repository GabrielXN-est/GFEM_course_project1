#include "lin_alg.h"
#include <thread>
#include <iostream>
#include <cmath>

bool verify{false};

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

void diagonalize(Matrix& K, Vector& F, std::size_t n)
{
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
}

// resolver sistema linear K * U = F
// método de eliminação de Gauss
Vector Gauss_elimination(Matrix K, Vector F)
{
    std::size_t n {F.size()};
    Vector U(n);

    if (verify)
    {
        std::thread t (verify_if_singular, K);
        diagonalize(K, F, n);
        t.join();
    }
    else
        diagonalize(K, F, n);

    for (std::size_t i {0}; i < n; i++)
    {
        U[i] = F[i]/K[i][i];
    }
    return std::move(U);
}

Matrix I (std::size_t size)
{
    Matrix identity(size, size);
    for (std::size_t i {0}; i < size; i++)
    {
        identity[i][i] = 1.0;
    }
    return identity;
}

double stop_condition(Vector&U, Vector&e, Matrix& K)
    {return (e.T()*K*e).determinant()/(U.T()*K*U).determinant();}

void solve_dependent_system(Matrix& K_, Vector& F_, double tol = std::pow(10, -12)) // Babuska et al.
{
    Matrix K (K_.mat.size(), K_[0].size());
    Matrix T (K_.mat.size(), K_[0].size());

    Vector F (F_.size()), U (F_.size()), e (F_.size()), r (F_.size());

    for (std::size_t i {0}; i < K_.mat.size(); i++)
    {
        for (std::size_t j {0}; j < K_[0].size(); j++)
        {
            K[i][j] = K_[i][j]/std::sqrt(K_[i][i] * K_[j][j]);
            if (i == j)
                T[i][j] = 1./std::sqrt(K_[i][i]);
            else
                T[i][j] = 0.;
        }
        F[i] = F_[i];
    }
 
    Matrix Ke {K + I(K.mat.size())};
    LU_factorization Ke_LU(Ke);

    int n_iter {0};
    do
    {
        Ke_LU.solve(F, U);
        r = F - K*U;
        Ke_LU.solve(r, e);
        
        n_iter++;
        if (n_iter > 1000)
            throw std::runtime_error("Warning: Maximum number of iterations reached without convergence. (" + std::to_string(n_iter) + ")");
    }
    while (stop_condition(U, e, K)> tol)
}

void LU_factorization::solve (Vector& F, Vector& U)
{
    std::size_t n {F.size()};
    Vector y(n);

    // Forward substitution Ly = F
    for (std::size_t i {0}; i < n; i++)
    {
        y[i] = F[i];
        for (std::size_t j {0}; j < i; j++)
        {
            y[i] -= LU[i][j] * y[j];
        }
    }

    // Backward substitution Ux = y
    for (std::size_t i {n-1}; i>= 0; i--)
    {
        U[i] = y[i];
        for (std::size_t j {i+1}; j < n; j++)
            {U[i] -= LU[i][j] * U[j];}
        U[i] /= LU[i][i];
    }
}