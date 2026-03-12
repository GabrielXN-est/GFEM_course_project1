#include "bars.h"

void lagrangian_bar::get_conectivity()
{
    conectivity = std::vector<int>(get_n_dofs ());

    for (std::size_t i {0}; i < conectivity.size(); i++)
        {conectivity[i] = Nod_list[i]->dofs[0];}
}

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
        if (node->dofs.size() == 0)
            node->dofs.push_back(dof0++);
    }
}

void p_hier_bar::get_conectivity()
{
    conectivity = std::vector<int>(get_n_dofs ());

    std::size_t dof {0}; 
    for (std::size_t i {0}; i < Nod_list.size(); i++)
    {
        for (std::size_t j {0};j< Nod_list[i]->dofs.size();j++)
            {conectivity[dof++] = Nod_list[i]->dofs[j];}
    }
}

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
        if (Nod_list[i]->dofs.size() == 0)
            if (i == 1)
                for (int j {0}; j < shape_func_order-1; j++)
                    {Nod_list[i]->dofs.push_back(dof0++);}
            else
                Nod_list[i]->dofs.push_back(dof0++);
    }
}