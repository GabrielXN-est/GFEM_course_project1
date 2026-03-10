#ifndef CONSTRAINED_BAR_H
#define CONSTRAINED_BAR_H
#include "lin_alg.h"
#include "node.h"
#include "Integration_points/integration_points.h"
#include "shape_functions.h"
#include "Inputs/read_input.h"

class Constrained_bar
{
    public:

    int id;
    int prop_id;
    double A {};
    std::vector<double> E {};
    std::vector<double> Exlim {};
    double C {};
    double L {};
    int shape_func_order {};

    Matrix K_local;
    Vector F_local;

    std::vector<int> N_list {};
    std::vector<Node*> Nod_list {};
    std::vector<int> conectivity {};

    body_functions* bf_func;

    Constrained_bar (int index, std::vector<int> nodeL, int prop) : 
        id{index}, prop_id{prop}, shape_func_order {static_cast<int>(nodeL.size()-1)},
        N_list {nodeL}, Nod_list (N_list.size()), conectivity (N_list.size()), 
        K_local (shape_func_order+1, shape_func_order+1), F_local (shape_func_order+1)
        {}    

    ~Constrained_bar() {delete bf_func;}
    
    // inicializa nós da barra
    void get_nodes(std::vector<Node>& nodevec)
    {
        int j {0};
        for (Node& i: nodevec)
        {
            if (std::count(N_list.begin(), N_list.end(), i.id))
                Nod_list[j++] = &i;
        }
        sortNodesByX(Nod_list);
    }

    // cria matriz de conectividade da barra
    void get_conectivity()
    {
        for (std::size_t i {0}; i < conectivity.size(); i++)
        {
            conectivity[i] = Nod_list[i]->dofs[0];
        }
    }

    // inicializa as propriedades da barra
    void get_properties(Properties prop)
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
        prop.bf_func = nullptr; // para evitar que o destrutor de prop delete o ponteiro de bf_func, já que a barra vai usar esse ponteiro
    }

    // do físico paro o mestre
    double mapping(double x, double xi, double Li)
    {
        return (x-xi)*2/Li;
    }
    // do mestre paro o físico
    double mapping_inv(double x, double xi, double Li)
    {
        return (x+1)*Li/2 + xi;
    }

    // integra tensões num intervalo da barra
    void integrate_B1_to_K(double Ei, double Li, double xi)
    {
        // para o trecho AE dNdx dNdxt
        integration_points ip1 {Gauss_quad_points(2*shape_func_order-2)};

        MDshape_functions_lag dNdx {shape_func_order};

        for (std::size_t pt_id {0}; pt_id < ip1.points.size(); pt_id++)
        {
            // obter dNdx e dNdxt
            dNdx(ip1.points[pt_id]);
            dNdx.mont_vector();

            // montar matriz de rigidez local
            K_local += (dNdx*dNdx.T())*(Li/2*ip1.weights[pt_id]*A*Ei);
        }
    }

    // integra tempo da mola distribuída
    void integrate_B2_to_K()
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

    // integra a força de corpo
    void integrate_BF_to_F()
    {
        integration_points ip {Gauss_quad_points(shape_func_order+bf_func->grau)};
        Mshape_functions_lag N {shape_func_order};

        for (std::size_t i {0}; i < ip.points.size(); i++)
        {
            N(ip.points[i]);
            N.mont_vector();
            
            F_local += N * static_cast<double>((*bf_func)(mapping_inv(ip.points[i], Nod_list[0]->x, L))*ip.weights[i]*L/2.);
        }
    }

    // calcula a matriz de rigidez local da barra
    void get_K()
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

    void start_el()
    {
        get_conectivity();
        get_K();
        integrate_BF_to_F();
    }
};

#endif