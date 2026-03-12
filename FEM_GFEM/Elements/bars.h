#ifndef CONSTRAINED_BAR_H
#define CONSTRAINED_BAR_H

#include "lin_alg.h"
#include "node.h"
#include "Integration_points/integration_points.h"
#include "shape_functions.h"
#include "Inputs/read_input.h"
#include "element.h"

class lagrangian_bar : public Bar
{
    public:
    lagrangian_bar (int index, std::vector<int> nodeL, int prop, int sh_o) : Bar {index, nodeL, prop, sh_o} {}   
    
    // getters
    // cria matriz de conectividade da barra
    void get_conectivity();

    shape_functions* get_shape_func();
    shape_functions* get_D_shape_func();

    void assign_dofs(int& dof0);
};

class p_hier_bar : public Bar
{
    public:
    p_hier_bar (int index, std::vector<int> nodeL, int prop, int sh_o) : Bar {index, nodeL, prop, sh_o} {}   
    
    // getters
    // cria matriz de conectividade da barra
    void get_conectivity();

    shape_functions* get_shape_func();
    shape_functions* get_D_shape_func();

    void assign_dofs(int& dof0);
};

#endif