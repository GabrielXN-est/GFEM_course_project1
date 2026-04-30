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
        dof = bc.no->dofs[static_cast<std::size_t>(bc.dof)];
        F_global[dof] = bc.value; 
    }

    K_global_pos = K_global;
    F_global_pos = F_global;

    for (BC_displacement& bc: bc_ds)
    {
        dof = bc.no->dofs[static_cast<std::size_t>(bc.dof)];
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
        idof = bc.no->dofs[static_cast<std::size_t>(bc.dof)];
        F_global[idof] = bc.value; 
    }

    std::vector<int> dof (bc_ds.size()); // dofs com bc de dirichilet

    for (std::size_t i {0}; i < bc_ds.size(); i++)
    {dof[i] = bc_ds[i].no->dofs[static_cast<std::size_t>(bc_ds[i].dof)];}

    K_global_pos = K_global;
    F_global_pos = F_global;

    K_global.clear();
    F_global.clear();

    int lin {0}, col {0};
    K_global = Matrix(K_global_pos.mat.size() - dof.size(), K_global_pos.mat.size() - dof.size());
    F_global = Vector(F_global_pos.size() - dof.size());
    T = Matrix (F_global.size(), F_global.size());

    for (std::size_t i {0}; i < K_global_pos.mat.size(); i++)
    {
        if (!(std::find(dof.begin(), dof.end(), static_cast<int>(i)) != dof.end()))
        {
            F_global[lin] = F_global_pos[i];
            for (BC_displacement  bc_dof : bc_ds)
            {
                F_global[lin] -= K_global_pos[i][bc_dof.no->dofs[static_cast<std::size_t>(bc_dof.dof)]] * bc_dof.value;
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

    for (std::size_t i{0}; i< T.mat.size(); i++)
        {T[i][i] = 1.;}
}

void Mesh::complete_U()
{
    std::vector<int> dof (bc_ds.size()); // dofs com bc de dirichilet
    std::vector<double> bcval (bc_ds.size()); // valores das bc de dirichilet

    for (std::size_t i {0}; i < bc_ds.size(); i++)
    {
        dof[i] = bc_ds[i].no->dofs[static_cast<std::size_t>(bc_ds[i].dof)];
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

void Mesh::create_scaled_global_system(bool get_condition_number)
{
    Matrix K_old {K_global};

    for (std::size_t i {0}; i < K_global.mat.size(); i++)
    {
        for (std::size_t j {0}; j < K_global[0].size(); j++)
        {
            K_global[i][j] = K_old[i][j]/std::sqrt(K_old[i][i] * K_old[j][j]);
            if (i == j)
                T[i][j] = 1./std::sqrt(K_old[i][i]);
            else
                T[i][j] = 0.;
        }
    }

    if (get_condition_number)
        scaled_condition_number = K_global.condition_number();

    F_global = T * F_global;
}

void Mesh::solve_dependent_system(double tol, int max_iter) // Babuska et al.
{
    Vector u (F_global.size()), e (F_global.size()), r (F_global.size());

    Matrix Ke {K_global + I(K_global.mat.size())};
    LU_factorization Ke_LU(Ke);

    Ke_LU.solve(F_global, u);

    int n_iter {0};
    do
    {
        u += e;
        r = F_global - K_global*u;
        Ke_LU.solve(r, e);
        
        n_iter++;

        if (n_iter > max_iter)
            throw std::runtime_error("Warning: Maximum number of iterations reached without convergence. (" + std::to_string(n_iter) + ") -- error of 10^-" + std::to_string(std::log10(stop_condition(u, e, K_global))));
    } while (stop_condition(u, e, K_global)> tol);

    U = T * u;
}

void Mesh::create_local_problem(double x0, double L, double nelem, std::string el_type, int pord, int E_pord, std::vector<double> geom_enr_r)
{ 
    //defaults
    if(!(el_type==""))
        el_type = eltype;
    if(pord == -1)
        pord = porder;
    if(E_pord == -1)
        E_pord = porder_Enrichment;
    if(geom_enr_r.size() == 1)
        geom_enr_r = geom_enr;

    if (el_type == "pGFEMBar_WD_S" || el_type == "pGFEMBar_2Proj")
        throw std::invalid_argument("Local problems not implemented for enrichments with non zero values on the nodes");

    // boundary conditions
    std::vector<double> new_f_bcs {}, new_d_bcs {};
    std::vector<int> new_f_bcs_pos {}, new_f_bcs_dofs {}, new_d_bcs_pos {}, new_d_bcs_dofs {};
    bool flag0 {false}, flagL {false};

    // Neuman boundary conditions in the borders of the local problem
    for (std::size_t i {0}; i < f_bcs_pos.size(); i++)
    {
        if ((f_bcs_pos[i] == 0) && (x0 == 0))
        {
            new_f_bcs.push_back(f_bcs[i]);
            new_f_bcs_pos.push_back(0);
            new_f_bcs_dofs.push_back(f_bcs_dofs[i]);
            flag0 = true;
        }
        else if ((f_bcs_pos[i] == 1) && (x0 + L == Lm))
        {
            new_f_bcs.push_back(f_bcs[i]);
            new_f_bcs_pos.push_back(1);
            new_f_bcs_dofs.push_back(f_bcs_dofs[i]);
            flagL = true;
        }
    }

    // dirichilet boundary conditions in the borders of the local problem
    for (std::size_t i {0}; i < d_bcs_pos.size(); i++)
    {
        if ((d_bcs_pos[i] == 0) && (x0 == 0))
        {
            new_d_bcs.push_back(d_bcs[i]);
            new_d_bcs_pos.push_back(0);
            new_d_bcs_dofs.push_back(d_bcs_dofs[i]);
            flag0 = true;
        }
        else if ((d_bcs_pos[i] == 1) && (x0 + L == Lm))
        {
            new_d_bcs.push_back(d_bcs[i]);
            new_d_bcs_pos.push_back(1);
            new_d_bcs_dofs.push_back(d_bcs_dofs[i]);
            flagL = true;
        }
    }

    // approximate boundary conditions if not in the border of theoriginal problem
    if (!flag0)
    {
        new_d_bcs.push_back(interpolate_solution(x0));
        new_d_bcs_pos.push_back(0);
        new_d_bcs_dofs.push_back(1);
    }
    if (!flagL)
    {
        new_d_bcs.push_back(interpolate_solution(x0+L));
        new_d_bcs_pos.push_back(1);
        new_d_bcs_dofs.push_back(1);
    }

    // name of the file to save
    std::string new_filealias {filealias + "_local_" + std::to_string(local_problems.size())};

    // create mesh and save in local problems vector
    local_problems.emplace_back(path, new_filealias, nelem, pord, el_type,
    L, E, Exlim, A, C,
    d_bcs, d_bcs_pos, d_bcs_dofs, // dirichilet boundary conditions
    f_bcs, f_bcs_pos, f_bcs_dofs, // Neumann Boundary conditions
    bf_func_id, alpha, xb, xi, 
    xgamma, E_pord, geom_enr_r);
}

double Mesh::interpolate_solution(double x)
{
    double u {0};
    for (Element* el : c_bars)
    {
        if (el->Nod_list[0]->x <= x && el->Nod_list[0]->x + el->el_size >= x)
        {
            shape_functions* sf {el->get_shape_func()};
            sf->operator()(el->mapping(x, el->Nod_list[0]->x, el->el_size));
            sf->mont_vector();
            Vector N (el->Ndof);
            el->Mont_N(N, sf, el->Nod_list, x, el->Ndof);
            for (std::size_t i {0}; i < N.vec.size(); i++)
                u += N.vec[i]*U[el->conectivity[i]];
            return u;
        }
    }
    throw std::out_of_range("Point x is out of the domain of the problem.");
    return u;
}

double Mesh::interpolate_D_of_solution(double x)
{
    double dudx {0};
    for (Element* el : c_bars)
    {
        if (el->Nod_list[0]->x <= x && el->Nod_list[0]->x + el->el_size >= x)
        {
            shape_functions* dsfdxiPoU {el->get_D_shape_func()};
            dsfdxiPoU->operator()(el->mapping(x, el->Nod_list[0]->x, el->el_size));
            dsfdxiPoU->mont_vector();
            Vector dNdx (el->Ndof);
            el->Mont_dNdx(dNdx, dsfdxiPoU, el->Nod_list, x, 2/el->el_size);

            for (std::size_t i {0}; i < dNdx.vec.size(); i++)
                dudx += dNdx.vec[i]*U[el->conectivity[i]];
            return dudx;
        }
    }
    throw std::out_of_range("Point x is out of the domain of the problem.");
    return dudx;
}

void Mesh::run()
{
    create_mesh();
    read_input(path+"/input_files/"+filealias+".txt", *this);
    assemble_direct();
    create_scaled_global_system(true);
    if (eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M" || eltype == "pGFEMBar" || eltype == "pGFEMBar_sc" || eltype == "pGFEMBar_2Proj" || eltype == "pSGFEMBar_2Proj")
        solve_dependent_system(std::pow(10, -30), 100000000);
    else
        solve();
    complete_U();
    local_problems.clear(); // tem de ser recomputados
}

// run local problems
void Mesh::run_local_problems()
{
    for (Mesh& local_mesh : local_problems)
    {
        local_mesh.run();
    }
}