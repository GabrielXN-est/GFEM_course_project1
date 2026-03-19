#include "bars.h"

// funções auxiliares
void assign_enrichment_dofs(Node* node, int& dof0)
{   for (size_t i {0}; i < node->enr.size(); i++)
        {node->dofs.push_back(dof0++);}}

// lagrangian_bar
shape_functions* lagrangian_bar::get_shape_func()
{
    return new Mshape_functions_lag(shape_func_order);
}

shape_functions* lagrangian_bar::get_D_shape_func()
{
    return new MDshape_functions_lag(shape_func_order);
}

void lagrangian_bar::assign_dofs(int& dof0)
{
    for (Node* node : Nod_list)
    {
        // PoU dofs
        if (node->dofs.size() == 0)
            node->dofs.push_back(dof0++);
        // Enrichment dofs
        assign_enrichment_dofs(node, dof0);
    }
}

// p_hier_bar
shape_functions* p_hier_bar::get_shape_func()
{
    return new Mshape_functions_p_hier(shape_func_order);
}

shape_functions* p_hier_bar::get_D_shape_func()
{
    return new MDshape_functions_p_hier(shape_func_order);
}

void p_hier_bar::assign_dofs(int& dof0)
{
    for (std::size_t i {0}; i < 3; i++)
    {
        // PoU dofs
        if (Nod_list[i]->dofs.size() == 0)
            if (i == 1)
                for (int j {0}; j < shape_func_order-1; j++)
                    {Nod_list[i]->dofs.push_back(dof0++);}
        // Enrichment dofs
        assign_enrichment_dofs(Nod_list[i], dof0);
    }
}

// p_GFEM_bar
void p_GFEM_bar::set_enrichments()
{
    for (Node* node : Nod_list)
    {
        if (node->p_enriched < Enrich)
        {
            for (int grau; grau < Enrich; grau++)
            {
                node->enr.push_back(new polinomial_enrichment_1D(-1, grau+1));;
            }
            node-> p_enriched = Enrich;
        }
    }
}

void p_GFEM_bar::start_el(std::vector<Node>& node_vec, int& dof0, std::vector<Properties>& pr_vec)
{
    get_nodes(node_vec);
    assign_dofs(dof0);
    set_enrichments();

    for (Properties& pr: pr_vec)
    {
        if (pr.id == prop_id)
        {
            get_properties(pr);
            break;
        }
    }

    // calcular matrizes locais
        get_conectivity();
        get_K();
        integrate_BF_to_F();
}