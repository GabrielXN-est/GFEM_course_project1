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
    lagrangian_bar (int index, std::vector<int> nodeL, int prop, int sh_o) : Bar {index, nodeL, prop, sh_o+1, sh_o} {}   
    
    // getters
    shape_functions* get_shape_func();
    shape_functions* get_D_shape_func();

    void assign_dofs(int& dof0);
};

class p_hier_bar : public Bar
{
    public:
    p_hier_bar (int index, std::vector<int> nodeL, int prop, int sh_o) : Bar {index, nodeL, prop, sh_o+1, sh_o} {}   
    
    // getters
    shape_functions* get_shape_func();
    shape_functions* get_D_shape_func();

    void assign_dofs(int& dof0);
};

class p_GFEM_bar : public lagrangian_bar
{
    public:
    int Enrich {};
    
    p_GFEM_bar (int index, std::vector<int> nodeL, int prop, int PoU_sh_o, int E_sh_o) :
     lagrangian_bar {index, nodeL, prop, (PoU_sh_o+1)*(E_sh_o+1)}, Enrich {E_sh_o} {}
    
    // definie enriquecimentos nos nós
    void set_enrichments();

    virtual void start_el(std::vector<Node>& node_vec, int& dof0, std::vector<Properties>& pr_vec) override;
};

#endif