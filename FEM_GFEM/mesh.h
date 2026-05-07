#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include "node.h"
#include "element.h"
#include "Bondeary_conditions.h"
#include "create_input.h"
#include <cmath>

// malha
class Mesh
{
    public:
    double Big_number {std::pow(10, 100)};
    double scaled_condition_number {};

    int nlocal_problems {0};
    std::vector<Mesh> local_problems {};

    std::vector<Node> nodes {};
    std::vector<Element*> c_bars {};
    std::vector<BC_displacement> bc_ds {};
    std::vector<BC_load> bc_l {};
    std::string name {};

    // criar essas matrizes
    Matrix K_global;
    Matrix K_global_pos;

    Matrix T; // scalling vector

    Vector F_global; 
    Vector F_global_pos; 

    Vector U;

    // Variaveis para gerar malha
    std::string path {};
    std::string filealias {};
    int nel {};
    int porder {};
    std::string eltype {};
    double Lm {};
    const std::vector<double>& E {};
    const std::vector<double>& Exlim {};
    double A {};
    double C {};
    const std::vector<double>& d_bcs {};
    const std::vector<int>& d_bcs_pos {};
    const std::vector<int>& d_bcs_dofs {};
    const std::vector<double>& f_bcs {};
    const std::vector<int>& f_bcs_pos {};
    const std::vector<int>& f_bcs_dofs {};
    int bf_func_id {};
    double alpha {};
    double xb {};
    double xi {};
    double xgamma {};
    int porder_Enrichment {};
    std::vector<double> geom_enr {};

    // constructors
    Mesh () {}

    Mesh (std::string Path, std::string filename, int nele, int pord, std::string elT,
    double l, const std::vector<double>& El, const std::vector<double>& Exliml, double a, double c,
    const std::vector<double>& D_bcs, const std::vector<int>& D_bcs_pos, const std::vector<int>& D_bcs_dofs, // dirichilet boundary conditions
    const std::vector<double>& F_bcs, const std::vector<int>& F_bcs_pos, const std::vector<int>& F_bcs_dofs, // Neumann Boundary conditions
    int bf_func, double alp=0.0, double x_b=0.0, double x_i = 0.0, 
    double xg = 0.0, int porder_E= 0, std::vector<double> geom_en = {}, int Nlocal_problems = 0) :
    path{Path}, filealias{filename}, nel{nele}, porder{pord}, eltype{elT},
    Lm{l}, E{El}, Exlim{Exliml}, A{a}, C{c},
    d_bcs{D_bcs}, d_bcs_pos{D_bcs_pos}, d_bcs_dofs{D_bcs_dofs},
    f_bcs{F_bcs}, f_bcs_pos{F_bcs_pos}, f_bcs_dofs{F_bcs_dofs},
    bf_func_id{bf_func}, alpha{alp}, xb{x_b}, xi{x_i}, 
    xgamma{xg}, porder_Enrichment{porder_E}, geom_enr{geom_en}, nlocal_problems{Nlocal_problems}
    {local_problems.reserve(Nlocal_problems);}

    ~Mesh () 
    {
        for (Element* el : c_bars)
            delete el;
    }

    //Setters
    void set_dofs(int dofs);

    void assign_nodes_biggest_vicinal_element_size();

    // Generate input file
    void create_mesh() 
    {
        generate_input(path+"/input_files/"+filealias+".txt", nel, porder, eltype,Lm,E, Exlim,A,C,
        d_bcs, d_bcs_pos, d_bcs_dofs, // dirichilet boundary conditions
        f_bcs, f_bcs_pos, f_bcs_dofs, // Neumann Boundary conditions
        bf_func_id, alpha, xb, xi, xgamma, porder_Enrichment, geom_enr, nlocal_problems);
    }

    // Asseblagem
    // BC pelo método da penalidade
    void assemble_penalty();

    // BC pelo método direto
    void assemble_direct();

    // Scale the stiffnes matrix
    void create_scaled_global_system(bool get_condition_number = false);

    // solve the system of equations
    void solve() 
    {
        U = T * Gauss_elimination(K_global, F_global);
    }
    void solve_dependent_system(double tol = std::pow(10, -12), int max_iter = 1000); // Babuska et al.

    // run the analysis
    void first_run();

    // run the analysis
    void run();
    
    // função para completar U com as condições de contorno se usado o método direto
    void complete_U ();

    double strain_energy();

    // Intepolate solution in a point x
    double interpolate_solution(double x);
    double interpolate_D_of_solution(double x);

    // cria um problema local ao repartir o domínio 
    void create_local_problem(double x0, double L, double nelem, std::string el_type="", int pord=-1, int E_pord=-1, std::vector<double> geom_enr_r = {1});

    // clear local problems to recompute them
    void clear_local_problems() {local_problems.clear();}

    // run local problems
    void run_local_problems();

    // reset the mesh elements (to run it again with new local problem enrichments)
    void reset_mesh_elements();
};

#endif