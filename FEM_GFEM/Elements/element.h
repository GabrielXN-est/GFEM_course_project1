#ifndef ELEMENT_H 
#define ELEMENT_H 

#include "lin_alg.h"
#include "node.h"
#include "integration_points.h"
#include "shape_functions.h"
#include "read_input.h"
#include "body_func.h"
#include "Enrichment.h"

class Element
{
    public:

    int id;
    int prop_id;
    int Ndof {}; // número de dofs
    double el_size {};
    Matrix K_local;
    Vector F_local;

    std::vector<int> N_list {};
    std::vector<Node*> Nod_list {};
    std::vector<int> conectivity {};

    Element (int index, std::vector<int> nodeL, int prop, int dofs) : 
        id{index},Ndof {dofs},
        N_list {nodeL}, Nod_list (N_list.size()), prop_id{prop},
        K_local (Ndof, Ndof), F_local (Ndof)
        {}    
    
    // getters
    // inicializa nós
    void get_nodes(std::vector<Node>& nodevec);

    // inicializa matrizes e vetores locais
    void Set_ndof();

        // inicliaza o elemento para assenblagem
        virtual void start_el(std::vector<Node>& node_vec, int& dof0, std::vector<Properties>& pr_vec)=0;
        virtual void start_local() = 0;
        // setters
            virtual void get_conectivity()=0;
            virtual void get_properties(Properties& prop)=0;
        // getters
            virtual void get_K()=0;
        // Mapeamento
            virtual double mapping(double x, double xi, double Li)=0;
            virtual double mapping_inv(double x, double xi, double Li)=0;
        //shape functions
            virtual shape_functions* get_shape_func()=0;
            virtual shape_functions* get_D_shape_func()=0;
            virtual void Mont_N(Vector& N, shape_functions* N_PoU, std::vector<Node*>& Nod_list, double x_real, int Ndof) = 0;
        // inicliaza os graus de liberdade dos nós do elemento
            virtual void assign_dofs(int& dof0)=0;
            
};

class Bar: public Element
{    
    public:
    int shape_func_order {}; // o grau da resposta do FEM puro
    int E_shape_func_order {0}; // o grau do maior enriquecimento polinomial

    double A {};
    std::vector<double> E {};
    std::vector<double> Exlim {};
    double C {};
    double L {};

    body_functions* bf_func;

    Bar (int index, std::vector<int> nodeL, int prop, int dofs, int sh_o, int hiper_dofs) : shape_func_order {sh_o}, Element (index, nodeL, prop, hiper_dofs){}    
    Bar (int index, std::vector<int> nodeL, int prop, int dofs, int sh_o) : shape_func_order {sh_o}, Element (index, nodeL, prop, dofs){}  
    ~Bar() {delete bf_func;}

    // inicializa as propriedades da barra
    void get_properties(Properties& prop);

    // determina o maior grau das funções de enriquecimento para determinação dos pontos de integração
    void assign_E_degree();

    // cria matriz de conectividade da barra
    void get_conectivity();

    //Integrais
    // integra tensões num intervalo da barra
    void integrate_B1_to_K(double Ei, double Li, double xi);

    // Cria o vetor das funções de forma do GFEM
    void Mont_N(Vector& N, shape_functions* N_PoU, std::vector<Node*>& Nod_list, double x_real, int Ndof);
    
    // integra termo da mola distribuída
    void integrate_B2_to_K();

    // integra a força de corpo
    void integrate_BF_to_F();

    // funções de mapeamento
    // do físico paro o mestre
    double mapping(double x, double xi, double Li){return (x-xi)*2/Li-1;}

    // do mestre paro o físico
    double mapping_inv(double x, double xi, double Li){return (x+1)*Li/2 + xi;}

    // calcula a matriz de rigidez local da barra
    void get_K();

    // inicliaza o elemento para assenblagem
    virtual void start_el(std::vector<Node>& node_vec, int& dof0, std::vector<Properties>& pr_vec);
    void start_local();
};

#endif