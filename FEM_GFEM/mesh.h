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

    void set_dofs(int dofs)
    {
        K_global = Matrix(dofs, dofs);
        F_global = Vector(dofs);
        U = Vector(dofs);
    }

    void assemble()
    {
        for (Constrained_bar& Cbar : c_bars)
        {
            for (std::size_t i {0}; i < Cbar.shape_func_order+1; i++)
            {
                F_global[Cbar.conectivity[i]] += Cbar.F_local[i];

                for (std::size_t j {0}; j < Cbar.shape_func_order+1; j++)
                {
                    K_global[Cbar.conectivity[i]][Cbar.conectivity[j]] += Cbar.K_local[i][j];
                }
            }
        }

        K_global_pos = K_global;

        int dof {};
        for (BC_displacement& bc: bc_ds)
        {
            dof = bc.no[0]->dofs[static_cast<std::size_t>(bc.dof)];
            K_global[dof][dof] = Big_number;
            F_global[dof] = Big_number * bc.value; 
        }

        // considerando cargas só nas laterais do elementos (cada ponto de aplicação de carga tem apenas 1 função de forma com valor não nulo)

        for (BC_load & bc: bc_l)
        {
            dof = bc.no[0]->dofs[static_cast<std::size_t>(bc.dof)];
            F_global[dof] = bc.value; 
        }
    }

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
        return 1/2 * (U.T() * K_global_pos * U).determinant();
    }
};

#endif