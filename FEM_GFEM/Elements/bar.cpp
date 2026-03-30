#include "element.h"

//Bar
void Bar::get_properties(Properties& prop)
{
    A = prop.A;
    C = prop.C;

    double xi {Nod_list[0]->x};
    L = el_size = Nod_list[N_list.size()-1]->x - xi;
    
    if (prop.E.size() == 1)
        E = prop.E;
    else
    {
        for (std::size_t i {0}; i < prop.E.size(); i++)
        {
            if (i==0 && prop.Exlim[i] > xi)
            {
                E.push_back(prop.E[i]);
                Exlim.push_back(prop.Exlim[i]);
            }
            else if (prop.Exlim[i-1] < xi && prop.Exlim[i] > xi ||
                     prop.Exlim[i-1] >= xi && prop.Exlim[i] <= xi + L ||
                     prop.Exlim[i-1] < xi + L && prop.Exlim[i] >= xi + L)
            {
                E.push_back(prop.E[i]);
                Exlim.push_back(prop.Exlim[i]);
            }
        }
    }
    bf_func = prop.bf_func->clone(); // para evitar que o destrutor de prop delete o ponteiro de bf_func, já que a barra vai usar esse ponteiro
}

void Bar::integrate_B1_to_K(double Ei, double Li, double xi)
{
    // para o trecho AE dNdx dNdxt
    integration_points ip1 {Gauss_quad_points(2*(shape_func_order+E_shape_func_order)-2)};

    shape_functions* dNdxiPoU {get_D_shape_func()};
    shape_functions* N_PoU {get_shape_func()};

    double dxidx {2/L}; // considerando mapeamento linear para a célula de integração
    double x_real {};

    for (std::size_t pt_id {0}; pt_id < ip1.points.size(); pt_id++)
    {
        x_real = mapping_inv(ip1.points[pt_id], xi, Li);
        // obter dNdx e dNdxt
            // avaliando o ponto de integração no sistema de coordenadas mestres
            // ao traduzir ele do sistema de cordenadas mestra da célula de integração para o físico e depois para o mestre do elemento
        dNdxiPoU->operator()(mapping(x_real, Nod_list[0]->x, L));
        dNdxiPoU->mont_vector();
        N_PoU->operator()(mapping(x_real, Nod_list[0]->x, L));
        N_PoU->mont_vector();

        // Consertar K
        int j {0};
        Vector dNdx (Ndof);
        for (std::size_t i {0}; i < Nod_list.size(); i++)
        {
            dNdx[j++] = dNdxiPoU->operator[](i) * dxidx;
            for (Enrichment* e : Nod_list[i]->enr)
            {
                dNdx[j++] = e->D(x_real)*N_PoU->operator[](i) + e->operator()(x_real)*dNdxiPoU->operator[](i)*dxidx;
            }
        }
        // montar matriz de rigidez local
        K_local += (dNdx*dNdx.T())*(Li/2*ip1.weights[pt_id]*A*Ei);
    }
}

void Bar::Mont_N(Vector& N, shape_functions* N_PoU, std::vector<Node*>& Nod_list, double x_real, int Ndof)
{
    int j {0};
    for (std::size_t i {0}; i < Nod_list.size(); i++)
    {
        N[j++] = N_PoU->operator[](i);
        for (Enrichment* e : Nod_list[i]->enr)
        {
            N[j++] = N_PoU->operator[](i)*e->operator()(x_real);
        }
    }
}

void Bar::integrate_B2_to_K()
{
    // para o trecho C N Nt
    integration_points ip2 {Gauss_quad_points(2*(shape_func_order+E_shape_func_order))};

    shape_functions* N_PoU {get_shape_func()};
    Vector N (Ndof);

    double x_real {};

    for (std::size_t pt_id {0}; pt_id < ip2.points.size(); pt_id++)
    {
        x_real = mapping_inv(ip2.points[pt_id], Nod_list[0]->x, L);
        
        // obter N
        N_PoU->operator()(ip2.points[pt_id]);
        N_PoU->mont_vector();
        Mont_N(N, N_PoU, Nod_list, x_real, Ndof);

        // montar matriz de rigidez local
        K_local += (N*N.T())*(L/2*ip2.weights[pt_id]*C);
    }
}

void Bar::integrate_BF_to_F()
{
    integration_points ip {Gauss_quad_points(shape_func_order+E_shape_func_order+bf_func->grau)};
    shape_functions* N_PoU {get_shape_func()};
    Vector N (Ndof);

    double x_real {};

    for (std::size_t i {0}; i < ip.points.size(); i++)
    {
        x_real = mapping_inv(ip.points[i], Nod_list[0]->x, L);

        N_PoU->operator()(ip.points[i]);
        N_PoU->mont_vector();
        Mont_N(N, N_PoU, Nod_list, x_real, Ndof);

        F_local += N * static_cast<double>((*bf_func)(x_real)*ip.weights[i]*L/2.);
    }
}

void Bar::get_K()
{
    integrate_B2_to_K();

    if (E.size() == 1)
        integrate_B1_to_K(E[0], L, Nod_list[0]->x);
    else
    {
        double xi {Nod_list[0]->x}, Li {0};
        double xi_ = xi;
        for (std::size_t i {0}; i < E.size(); i++)
        {
            if (Exlim[i] >= xi_)
            {
                Li = (Exlim[i] < xi_ + L) ? Exlim[i] - xi : L - xi + xi_;
                integrate_B1_to_K(E[i], Li, xi);
                xi = Exlim[i];
            }
            if (xi >= L + xi_)
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
}

void Bar::start_local()
{
    // calcular matrizes locais
    Set_ndof();
    assign_E_degree();
    get_conectivity();
    get_K();
    integrate_BF_to_F();
}

void Bar::assign_E_degree()
{
    for (Node* node : Nod_list)
    {
        for (Enrichment* e : node->enr)
        {
            if (e->grau > E_shape_func_order)
                E_shape_func_order = e->grau;
        }
    }
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