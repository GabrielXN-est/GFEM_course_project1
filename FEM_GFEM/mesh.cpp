#include "mesh.h"

void Mesh::set_dofs(int dofs)
{
    K_global = Matrix(dofs, dofs);
    F_global = Vector(dofs);
    U = Vector(dofs);
}

// Asseblagem
void Mesh::assemble_penalty()
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

//void Mesh::assemble_direct()
// implementar


