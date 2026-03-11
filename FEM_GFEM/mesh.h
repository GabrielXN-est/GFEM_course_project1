#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include "node.h"
#include "constrained_bar.h"
#include "Bondeary_conditions.h"
#include <cmath>

// malha
class Mesh
{
    public:
    double Big_number {std::pow(10, 100)};

    std::vector<Node> nodes {};
    std::vector<Constrained_bar> c_bars {};
    std::vector<BC_displacement> bc_ds {};
    std::vector<BC_load> bc_l {};
    std::string name {};

    // criar essas matrizes
    Matrix K_global;
    Matrix K_global_pos;
    Vector F_global; 

    Vector U;

    Mesh () {}

    //Setters
    void set_dofs(int dofs);

    // Asseblagem
    // BC pelo método da penalidade
    void assemble_penalty();

    // BC pelo método direto
    void assemble_direct();

    void solve()
    {
        Matrix K {K_global};
        Vector F {F_global};
        // resolver sistema linear K_global * U = F_global
        // método de eliminação de Gauss
        int n {static_cast<int>(F_global.size())};

        for (int i {0}; i < n; i++)
        {
            for (int j {0}; j < n; j++)
            {
                if (j != i)
                {
                    double factor {K_global[j][i]/K_global[i][i]};
                    for (int k {i}; k < n; k++)
                    {
                        K[j][k] -= factor * K_global[i][k];
                    }
                    F[j] -= factor * F_global[i];
                }
            }
        }

        for (int i {0}; i < n; i++)
        {
            U[i] = F[i]/K[i][i];
        }
    }

    double energy_norm()
    {
        Matrix temp1 {K_global_pos * U};
        return 1./2. * (U.T() * (K_global_pos * U)).determinant();
    }
};

#endif