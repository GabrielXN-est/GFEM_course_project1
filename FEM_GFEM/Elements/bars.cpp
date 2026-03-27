#include "bars.h"

// funções auxiliares
void assign_enrichment_dofs(Node* node, int& dof0)
{   
    int lacking_dofs {node->enr.size() + 1 - node->dofs.size()};
    for (int i {0}; i < lacking_dofs; i++)
        {node->dofs.push_back(dof0++);}
}

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
        if (node->dofs.size() < node->enr.size() + 1)
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
            for (int grau {0}; grau < Enrich; grau++)
            {
                node->enr.push_back(new polinomial_enrichment_1D(-1, grau+1, shifted, scaled, node));;
            }
            node-> p_enriched = Enrich;
        }
    }
}

void p_GFEM_bar::start_el(std::vector<Node>& node_vec, int& dof0, std::vector<Properties>& pr_vec)
{
    get_nodes(node_vec);
    set_enrichments();
    assign_dofs(dof0);

    for (Properties& pr: pr_vec)
    {
        if (pr.id == prop_id)
        {
            get_properties(pr);
            break;
        }
    }
}

// p_GFEM_bar_weak_disc
void p_GFEM_bar_weak_disc::set_enrichments_desc()
{
    set_enrichments();
    for (Node* node : Nod_list)
    {
        if (Enrich == 0 && node-> p_enriched >= 0)
            enrichment->assign_to_node(*node);
        else if (node->ep_enriched < Enrich)
        {
            enrichment->assign_to_node(*node);
            for (int grau {0}; grau < Enrich; grau++)
            {
                node->enr.push_back(new Pair_enrichment_1D(
                    new polinomial_enrichment_1D(-1, grau+1, shifted, scaled, node),
                    enrichment->create_copy(*node)
                    ));
            }
        }
        node-> ep_enriched = Enrich;
    }
}

void p_GFEM_bar_weak_disc::start_el(std::vector<Node>& node_vec, int& dof0, std::vector<Properties>& pr_vec)
{
    get_nodes(node_vec);
    set_enrichments_desc();
    assign_dofs(dof0);

    for (Properties& pr: pr_vec)
    {
        if (pr.id == prop_id)
        {
            get_properties(pr);
            break;
        }
    }
}