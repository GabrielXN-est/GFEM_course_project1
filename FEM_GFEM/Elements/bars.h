#ifndef CONSTRAINED_BAR_H
#define CONSTRAINED_BAR_H

#include "element.h"

class lagrangian_bar : public Bar
{
    public:
    
    lagrangian_bar (int index, std::vector<int> nodeL, int prop, int sh_o, int hiper_dofs) : Bar {index, nodeL, prop, sh_o+1, sh_o, hiper_dofs} {} 
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
     lagrangian_bar {index, nodeL, prop, PoU_sh_o, E_sh_o + PoU_sh_o + 2}, 
     Enrich {E_sh_o}
     {E_shape_func_order = Enrich;};

    p_GFEM_bar (int index, std::vector<int> nodeL, int prop, int PoU_sh_o, int E_sh_o, int hiper_dofs) : 
     lagrangian_bar {index, nodeL, prop, PoU_sh_o, hiper_dofs}, 
     Enrich {E_sh_o}
     {E_shape_func_order = Enrich;};
    
    // definie enriquecimentos nos nós
    void set_enrichments();

    virtual void start_el(std::vector<Node>& node_vec, int& dof0, std::vector<Properties>& pr_vec) override;
};

class p_GFEM_bar_weak_disc : public p_GFEM_bar
{
    public:
    Enrichment* enrichment {};
    
    p_GFEM_bar_weak_disc (int index, std::vector<int> nodeL, int prop, int PoU_sh_o, int E_sh_o, Enrichment* enr) : 
     p_GFEM_bar (index, nodeL, prop, PoU_sh_o, E_sh_o, 2*E_sh_o + PoU_sh_o + 2), enrichment {enr} 
     {E_shape_func_order = 2*Enrich+enr->grau;}

    ~p_GFEM_bar_weak_disc() {delete enrichment;}

    // definie enriquecimentos nos nós
    void set_enrichments_desc();
};

#endif