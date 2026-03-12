#ifndef ELEMENT_H 
#define ELEMENT_H 

#include "lin_alg.h"
#include "node.h"
#include "integration_points.h"
#include "shape_functions.h"
#include "read_input.h"
#include "body_func.h"

class Element
{
    public:

    int id;
    int prop_id;
    int shape_func_order {};

    Matrix K_local;
    Vector F_local;

    std::vector<int> N_list {};
    std::vector<Node*> Nod_list {};
    std::vector<int> conectivity {};

    Element (int index, std::vector<int> nodeL, int prop, int sh_o) : 
        id{index},shape_func_order {sh_o},
        N_list {nodeL}, Nod_list (N_list.size()), prop_id{prop},
        K_local (shape_func_order+1, shape_func_order+1), F_local (shape_func_order+1)
        {}    
    
    // getters
    // inicializa nós
    void get_nodes(std::vector<Node>& nodevec);

    // inicliaza os graus de liberdade dos nós do elemento
    virtual void assign_dofs(int& dof0)=0;

    //funções abstratas
    // cria matriz de conectividade 
    virtual void get_conectivity()=0;

    // inicializa as propriedades
    virtual void get_properties(Properties& prop)=0;

    // calcula a matriz de rigidez local da barra
    virtual void get_K()=0;

    // inicliaza o elemento para assenblagem
    virtual void start_el()=0;
};

class Bar: public Element
{
    public:

    double A {};
    std::vector<double> E {};
    std::vector<double> Exlim {};
    double C {};
    double L {};

    body_functions* bf_func;

    Bar (int index, std::vector<int> nodeL, int prop, int sh_o) : Element (index, nodeL, prop, sh_o){}    

    ~Bar() {delete bf_func;}
    
    // getters
    // cria matriz de conectividade da barra
    virtual void get_conectivity()=0;

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

    // Auxiliares
    int get_n_dofs();

    virtual shape_functions* get_shape_func()=0;
    virtual shape_functions* get_D_shape_func()=0;
};
#endif