#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include "node.h"
#include "element.h"
#include "Bondeary_conditions.h"
#include <cmath>

// malha
class Mesh
{
    public:
    double Big_number {std::pow(10, 100)};

    std::vector<Node> nodes {};
    std::vector<Element*> c_bars {};
    std::vector<BC_displacement> bc_ds {};
    std::vector<BC_load> bc_l {};
    std::string name {};

    // criar essas matrizes
    Matrix K_global;
    Matrix K_global_pos;
    Vector F_global; 
    Vector F_global_pos; 

    Vector U;

    Mesh () {}

    ~Mesh () 
    {
        for (Element* el : c_bars)
            delete el;
    }

    //Setters
    void set_dofs(int dofs);

    void assign_nodes_biggest_vicinal_element_size();

    // Asseblagem
    // BC pelo método da penalidade
    void assemble_penalty();

    // BC pelo método direto
    void assemble_direct();

    void solve() {U = Gauss_elimination(K_global, F_global);}
    void solve_dependent_system(double tol = std::pow(10, -12), int max_iter = 1000); // Babuska et al.

    // função para completar U com as condições de contorno se usado o método direto
    void complete_U ();

    double strain_energy();
};

#endif