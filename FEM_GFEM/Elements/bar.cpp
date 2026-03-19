#include "element.h"

//Bar
void Bar::get_properties(Properties& prop)
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

void Bar::integrate_B1_to_K(double Ei, double Li, double xi)
{
    // para o trecho AE dNdx dNdxt
    integration_points ip1 {Gauss_quad_points(2*(shape_func_order+E_shape_func_order)-2)};

    shape_functions* dNdxi {get_D_shape_func()};

    double dxidx {2/Li}; // considerando mapeamento linear para a célula de integração
    double x_real {};

    for (std::size_t pt_id {0}; pt_id < ip1.points.size(); pt_id++)
    {
        // obter dNdx e dNdxt
            // avaliando o ponto de integração no sistema de coordenadas mestres
            // ao traduzir ele do sistema de cordenadas mestra da célula de integração para o físico e depois para o mestre do elemento
        dNdxi->operator()(mapping(mapping_inv(ip1.points[pt_id], xi, Li), Nod_list[0]->x, L));
        
        // montar matriz de rigidez local
        K_local += ((*dNdxi)*(*dNdxi).T())*(dxidx*dxidx*Li/2*ip1.weights[pt_id]*A*Ei);
    }
}

void Bar::integrate_B2_to_K()
{
    // para o trecho C N Nt
    integration_points ip2 {Gauss_quad_points(2*(shape_func_order+E_shape_func_order))};

    shape_functions* N {get_shape_func()};

    for (std::size_t pt_id {0}; pt_id < ip2.points.size(); pt_id++)
    {
        // obter N
        N->operator()(ip2.points[pt_id]);
        N->mont_vector();

        // montar matriz de rigidez local
        K_local += ((*N)*(*N).T())*(L/2*ip2.weights[pt_id]*C);
    }
}

void Bar::integrate_BF_to_F()
{
    integration_points ip {Gauss_quad_points(shape_func_order+E_shape_func_order+bf_func->grau)};
    shape_functions* N {get_shape_func()};

    for (std::size_t i {0}; i < ip.points.size(); i++)
    {
        N->operator()(ip.points[i]);
        N->mont_vector();

        double x {(*bf_func)(mapping_inv(ip.points[i], Nod_list[0]->x, L))};

        F_local += (*N) * static_cast<double>((*bf_func)(mapping_inv(ip.points[i], Nod_list[0]->x, L))*ip.weights[i]*L/2.);
    }
}

void Bar::get_K()
{
    integrate_B2_to_K();

    if (E.size() == 1)
        integrate_B1_to_K(E[0], L, Nod_list[0]->x);
    else
    {
        double xi {Nod_list[0]->x};
        for (std::size_t i {0}; i < E.size(); i++)
        {
            integrate_B1_to_K(E[i], Exlim[i]-xi, xi);
            xi = Exlim[i];
            break;
        }
    }
}

void Bar::start_el(std::vector<Node>& node_vec, int& dof0, std::vector<Properties>& pr_vec)
{
    get_nodes(node_vec);
    assign_dofs(dof0);
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

void Bar::get_conectivity()  // dofs ordenados seguindo a ordem dos nós
{
    conectivity = std::vector<int>(Ndof);

    std::size_t dof {0}; 
    for (std::size_t i {0}; i < Nod_list.size(); i++)
    {
        for (std::size_t j {0};j< Nod_list[i]->dofs.size();j++)
            {conectivity[dof++] = Nod_list[i]->dofs[j];}
    }
}