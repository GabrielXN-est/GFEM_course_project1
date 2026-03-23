#ifndef ENRICHMENT_H
#define ENRICHMENT_H

#include "node.h"

class Enrichment
{
    public:
    int id {};
    int grau {};
    Enrichment(int index, int g): id{index}, grau {g}{};

    virtual void assign_to_node(Node& node);

    //funções abstratas
        virtual Enrichment* create_copy(Node& node) =0 ;
        // Avaliadores da função
        virtual double operator()(double x) = 0;
        virtual double D(double x) = 0;
};

class Sukumar_enrichment_1D : public Enrichment
{
    public:
    double xGamma {};
    int grau {1};

    Sukumar_enrichment_1D(int index, double xg) : Enrichment{index, 1}, xGamma {xg} {};

    double operator()(double x){return std::abs(x-xGamma);}

    double D(double x);

    Enrichment* create_copy(Node& node);
};

class Moes_enrichment_1D : public Enrichment
{
    public:
    double xGamma {};
    int grau {1};
    Node* node_ptr {nullptr};

    Moes_enrichment_1D(int index, double xg) : Enrichment{index, 1}, xGamma {xg} {};
    Moes_enrichment_1D(int index, double xg, Node* node) : Enrichment{index, 1}, xGamma {xg}, node_ptr {node} {};

    double operator()(double x);

    double D(double x);

    Enrichment* create_copy(Node& node);
};

class polinomial_enrichment_1D : public Enrichment
{
    public:
    double grau {};
    bool shifted {};
    bool scaled {}; // se scaled = true, executar nos nós a função Mesh::assign_nodes_biggest_vicinal_element_size();
    Node* node_ptr {nullptr};

    polinomial_enrichment_1D(int index, double g, bool sh = true, bool sc = false) :
        Enrichment{index, g}, grau {g}, shifted {sh}, scaled {sc} {};
    polinomial_enrichment_1D(int index, double g, bool sh, bool sc, Node* node) :
        Enrichment{index, g}, grau {g}, shifted {sh}, scaled {sc}, node_ptr {node}{};

    Enrichment* create_copy(Node& node);

    double operator()(double x)
    {
        double temp {x};
        if (shifted)
            temp -= node_ptr->x;
        if (scaled)            
            temp /= node_ptr->biggest_vicinal_element_size;
        return std::pow(temp, grau);
    }

    double D(double x)
    {
        double temp {x};
        if (shifted)
            temp -= node_ptr->x;
        if (scaled)   
        {    
            temp /= node_ptr->biggest_vicinal_element_size;
            return grau*std::pow(temp, grau-1)/ node_ptr->biggest_vicinal_element_size;
        }
        return grau*std::pow(temp, grau-1);
    }
};

// produto entre 2 enriquecimentos
class Pair_enrichment_1D : public Enrichment
{
    public:
    Enrichment* enr1 {nullptr};
    Enrichment* enr2 {nullptr};

    Pair_enrichment_1D(Enrichment* enr1, Enrichment* enr2) : enr1 {enr1}, enr2 {enr2},
    Enrichment{0, enr1->grau+enr2->grau} {};
    ~Pair_enrichment_1D() 
    {
        delete enr1;
        delete enr2;
    };

    double operator()(double x) {return (*enr1)(x) * (*enr2)(x);}

    double D(double x) {return enr1->D(x) * (*enr2)(x) + (*enr1)(x) * enr2->D(x);}

    Enrichment* create_copy(Node& node);
};

// Função para ordenar um vetor de nodes pela posição x
inline void sort_Enr_by_id(std::vector<Enrichment*>& nodes)
{
    std::sort(nodes.begin(), nodes.end(), 
              [](const Enrichment* a, const Enrichment  * b) {
                  return a->id < b->id;
              });
}
#endif // ENRICHMENT_H