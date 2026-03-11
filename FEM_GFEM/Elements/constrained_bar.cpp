#include "constrained_bar.h"

void Constrained_bar::get_nodes(std::vector<Node>& nodevec)
{
    int j {0};
    for (Node& i: nodevec)
    {
        if (std::count(N_list.begin(), N_list.end(), i.id))
            Nod_list[j++] = &i;
    }
    sortNodesByX(Nod_list);
}

void Constrained_bar::get_conectivity()
{
    for (std::size_t i {0}; i < conectivity.size(); i++)
    {
        conectivity[i] = Nod_list[i]->dofs[0];
    }
}

void Constrained_bar::get_properties(Properties& prop)
{
    A = prop.A;
    C = prop.C;

    double xi {Nod_list[0]->x};
    L = Nod_list[N_list.size()-1]->x - xi;
    
    if (prop.E.size() == 1)
        E = prop.E;
    else
    {
        for (std::size_t i {0}; i < prop.E.size(); i++)
        {
            if (prop.Exlim[i] > xi && prop.Exlim[i] <= xi + L)
            {
                E.push_back(prop.E[i]);
                Exlim.push_back(prop.Exlim[i]);
            }
        }
    }
    bf_func = prop.bf_func;
    prop.bf_func = prop.bf_func->clone(); // para evitar que o destrutor de prop delete o ponteiro de bf_func, já que a barra vai usar esse ponteiro
}

void Constrained_bar::integrate_B1_to_K(double Ei, double Li, double xi)
{
    // para o trecho AE dNdx dNdxt
    integration_points ip1 {Gauss_quad_points(2*shape_func_order-2)};

    MDshape_functions_lag dNdxi {shape_func_order};

    double dxidx {2/Li}; // considerando mapeamento linear

    for (std::size_t pt_id {0}; pt_id < ip1.points.size(); pt_id++)
    {
        // obter dNdx e dNdxt
        dNdxi(ip1.points[pt_id]);
        dNdxi.mont_vector();

        // montar matriz de rigidez local
        K_local += (dNdxi*dNdxi.T())*(dxidx*dxidx*Li/2*ip1.weights[pt_id]*A*Ei);
    }
}

void Constrained_bar::integrate_B2_to_K()
{
    // para o trecho C N Nt
    integration_points ip2 {Gauss_quad_points(2*shape_func_order)};

    Mshape_functions_lag N {shape_func_order};

    for (std::size_t pt_id {0}; pt_id < ip2.points.size(); pt_id++)
    {
        // obter N
        N(ip2.points[pt_id]);
        N.mont_vector();

        // montar matriz de rigidez local
        K_local += (N*N.T())*(L/2*ip2.weights[pt_id]*C);
    }
}

void Constrained_bar::integrate_BF_to_F()
{
    integration_points ip {Gauss_quad_points(shape_func_order+bf_func->grau)};
    Mshape_functions_lag N {shape_func_order};

    for (std::size_t i {0}; i < ip.points.size(); i++)
    {
        N(ip.points[i]);
        N.mont_vector();

        double x {(*bf_func)(mapping_inv(ip.points[i], Nod_list[0]->x, L))};

        F_local += N * static_cast<double>((*bf_func)(mapping_inv(ip.points[i], Nod_list[0]->x, L))*ip.weights[i]*L/2.);
    }
}

void Constrained_bar::get_K()
{
    integrate_B2_to_K();

    if (E.size() == 1)
    {
        integrate_B1_to_K(E[0], L, Nod_list[0]->x);
    }
    else
    {
        double temp {Nod_list[0]->x};
        for (std::size_t i {0}; i < E.size(); i++)
        {
            integrate_B1_to_K(E[i], Exlim[i]-temp, Nod_list[0]->x);
            temp = Exlim[i];
            break;
        }
    }
}

void Constrained_bar::start_el()
{
    get_conectivity();
    get_K();
    integrate_BF_to_F();
}
