#include "mesh.h"
#include <iostream>

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
        for (std::size_t i {0}; i < el->Ndof; i++)
        {
            F_global[el->conectivity[i]] += el->F_local[i];

            for (std::size_t j {0}; j < el->Ndof; j++)
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
        for (std::size_t i {0}; i < el->Ndof; i++)
        {
            F_global[el->conectivity[i]] += el->F_local[i];

            for (std::size_t j {0}; j < el->Ndof; j++)
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
            for (BC_displacement  bc_dof : bc_ds)
            {
                
                double temp {K_global_pos[i][bc_dof.dof] * bc_dof.value};
                F_global[lin] -= K_global_pos[i][bc_dof.dof] * bc_dof.value;
            }
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

double Mesh::strain_energy()
{
    return 1./2. * (U.T() * (K_global_pos * U)).determinant();
}

double stop_condition(Vector&U, Vector&e, Matrix& K)
    {return (e.T()*K*e).determinant()/(U.T()*K*U).determinant();}

void Mesh::solve_dependent_system(double tol, int max_iter) // Babuska et al.
{
    Matrix K (K_global.mat.size(), K_global[0].size());
    Matrix T (K_global.mat.size(), K_global[0].size());

    Vector F (F_global.size()), u (F_global.size()), e (F_global.size()), r (F_global.size());

    for (std::size_t i {0}; i < K_global.mat.size(); i++)
    {
        for (std::size_t j {0}; j < K_global[0].size(); j++)
        {
            K[i][j] = K_global[i][j]/std::sqrt(K_global[i][i] * K_global[j][j]);
            if (i == j)
                T[i][j] = 1./std::sqrt(K_global[i][i]);
            else
                T[i][j] = 0.;
        }
        F[i] = F_global[i];
    }

    F = T * F;
 
    Matrix Ke {K + I(K.mat.size())};
    LU_factorization Ke_LU(Ke);

    Ke_LU.solve(F, u);

    int n_iter {0};
    do
    {
        u += e;
        r = F - K*u;
        Ke_LU.solve(r, e);
        
        n_iter++;

        if (n_iter > max_iter)
            throw std::runtime_error("Warning: Maximum number of iterations reached without convergence. (" + std::to_string(n_iter) + ")");
    } while (stop_condition(u, e, K)> tol);

    U = T * u;
}