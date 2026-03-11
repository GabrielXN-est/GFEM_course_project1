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
    
    // getters
    // inicializa nós da barra
    void get_nodes(std::vector<Node>& nodevec);

    // cria matriz de conectividade da barra
    void get_conectivity();

    // inicializa as propriedades da barra
    void get_properties(Properties& prop);

    // funções de mapeamento
    // do físico paro o mestre
    double mapping(double x, double xi, double Li){return (x-xi)*2/Li;}

    // do mestre paro o físico
    double mapping_inv(double x, double xi, double Li){return (x+1)*Li/2 + xi;}

    // Cálculo de K e F
    // integra tensões num intervalo da barra
    void integrate_B1_to_K(double Ei, double Li, double xi);

    // integra termo da mola distribuída
    void integrate_B2_to_K();

    // integra a força de corpo
    void integrate_BF_to_F();

    // calcula a matriz de rigidez local da barra
    void get_K();

    // inicliaza o elemento para assenblagem
    void start_el();
};

#endif