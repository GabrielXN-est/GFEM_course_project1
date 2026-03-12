#include "mesh.h"

void Mesh::set_dofs(int dofs)
{
    K_global = Matrix(dofs, dofs);
    F_global = Vector(dofs);
}

// Asseblagem
void Mesh::assemble_penalty()
{
    for (Element* el : c_bars)
    {
        for (std::size_t i {0}; i < el->shape_func_order+1; i++)
        {
            F_global[el->conectivity[i]] += el->F_local[i];

            for (std::size_t j {0}; j < el->shape_func_order+1; j++)
            {
                K_global[el->conectivity[i]][el->conectivity[j]] += el->K_local[i][j];
            }
        }
    }
    
    // considerando cargas só nas laterais do elementos (cada ponto de aplicação de carga tem apenas 1 função de forma com valor não nulo)
    int dof {};
    for (BC_load & bc: bc_l)
    {
        dof = bc.no[0]->dofs[static_cast<std::size_t>(bc.dof)];
        F_global[dof] = bc.value; 
    }

    K_global_pos = K_global;
    F_global_pos = F_global;

    for (BC_displacement& bc: bc_ds)
    {
        dof = bc.no[0]->dofs[static_cast<std::size_t>(bc.dof)];
        K_global[dof][dof] = Big_number;
        F_global[dof] = Big_number * bc.value; 
    }
}

void Mesh::assemble_direct()
{
    for (Element* el : c_bars)
    {
        for (std::size_t i {0}; i < el->shape_func_order+1; i++)
        {
            F_global[el->conectivity[i]] += el->F_local[i];

            for (std::size_t j {0}; j < el->shape_func_order+1; j++)
            {
                K_global[el->conectivity[i]][el->conectivity[j]] += el->K_local[i][j];
            }
        }
    }

    // considerando cargas só nas laterais do elementos (cada ponto de aplicação de carga tem apenas 1 função de forma com valor não nulo)
    int idof {};
    for (BC_load & bc: bc_l)
    {
        idof = bc.no[0]->dofs[static_cast<std::size_t>(bc.dof)];
        F_global[idof] = bc.value; 
    }

    std::vector<int> dof (bc_ds.size()); // dofs com bc de dirichilet

    for (std::size_t i {0}; i < bc_ds.size(); i++)
    {dof[i] = bc_ds[i].no[0]->dofs[static_cast<std::size_t>(bc_ds[i].dof)];}

    K_global_pos = K_global;
    F_global_pos = F_global;

    K_global.clear();
    F_global.clear();

    int lin {0}, col {0};
    K_global = Matrix(K_global_pos.mat.size() - dof.size(), K_global_pos.mat.size() - dof.size());
    F_global = Vector(F_global_pos.size() - dof.size());

    for (std::size_t i {0}; i < K_global_pos.mat.size(); i++)
    {
        if (!(std::find(dof.begin(), dof.end(), static_cast<int>(i)) != dof.end()))
        {
            F_global[lin] = F_global_pos[i];
            for (std::size_t j {0}; j < K_global_pos.mat[i].size(); j++)
            {
                if (!(std::find(dof.begin(), dof.end(), static_cast<int>(j)) != dof.end()))
                {
                    K_global[lin][col++] = K_global_pos[i][j];
                }
            }
            col = 0;
            lin++;
        }
    }
}

void Mesh::complete_U()
{
    std::vector<int> dof (bc_ds.size()); // dofs com bc de dirichilet
    std::vector<double> bcval (bc_ds.size()); // valores das bc de dirichilet

    for (std::size_t i {0}; i < bc_ds.size(); i++)
    {
        dof[i] = bc_ds[i].no[0]->dofs[static_cast<std::size_t>(bc_ds[i].dof)];
        bcval[i] = bc_ds[i].value;
    }

    Vector U_temp (K_global_pos.mat.size());
    int lin {0};
    std::vector<int>::iterator it;

    for (std::size_t i {0}; i< K_global_pos.mat.size(); i++)
    {
        it = std::find(dof.begin(), dof.end(), static_cast<int>(i));
        if (it != dof.end())
        {
            U_temp[i] = bcval[std::distance(dof.begin(), it)];
        }
        else
        {
            U_temp[i] = U[lin++];
        }
    }
    U = U_temp;
}

double Mesh::energy_norm()
{
    Matrix temp1 {K_global_pos * U};
    return 1./2. * (U.T() * (K_global_pos * U)).determinant();
}