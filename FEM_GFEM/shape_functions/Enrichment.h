#ifndef ENRICHMENT_H
#define ENRICHMENT_H

#include "node.h"

class Enrichment
{
    public:
    int id {};
    int grau {};
    Enrichment(int index, int g): id{index}, grau {g}{};

    virtual void assign_to_node(Node& node) =0 ;
};

class Sukumar_enrichment_1D : public Enrichment
{
    public:
    double xGamma {};
    int grau {1};

    Sukumar_enrichment_1D(int index, double xg) : Enrichment{index, 1}, xGamma {xg} {};

    void assign_to_node(Node& node);
};

class polinomial_enrichment_1D : public Enrichment
{
    public:
    double grau {};

    polinomial_enrichment_1D(int index, double g) : Enrichment{index, g}, grau {g} {};

    void assign_to_node(Node& node);
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