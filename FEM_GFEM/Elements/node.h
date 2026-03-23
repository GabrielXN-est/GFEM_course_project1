#ifndef NODE_H
#define NODE_H

#include <vector>
#include <algorithm>
#include "Enrichment.h"

// agrupamento de dofs
class Node
{
    public:
    int id;
    double x {};
    std::vector<int> dofs {};
    std::vector<int> enr_ids {};
    std::vector<Enrichment*> enr {};
    std::vector<Element*> vicinal_elements {};

    int p_enriched {0}; // flag para indicar se o nó é enriquecido por polinômios até determinado grau
    int polinomial_order_of_enrichment {0}; // indica o grau do polinômio de enriquecimento do nó
    double biggest_vicinal_element_size {}; // para enriquecimentos polinomiais, indica o tamanho do maior elemento vicinal ao nó

    Node (int index, double x_coord, std::vector<int> enrichment_ids) : id{index}, x {x_coord}, enr_ids{enrichment_ids}{}
    Node ();
    
    ~Node ()
    {
        for (Enrichment* e : enr)
            delete e;
    }

    void get_polinomial_order_of_enrichment()
    {
        for (Enrichment* e : enr)
        {
            if (e -> grau > polinomial_order_of_enrichment)
                polinomial_order_of_enrichment = e -> grau;
        }
    }
};

// Função para ordenar um vetor de nodes pela posição x
inline void sortNodesByX(std::vector<Node*>& nodes)
{
    std::sort(nodes.begin(), nodes.end(), 
              [](const Node* a, const Node* b) {
                  return a->x < b->x;
              });
}

inline void sortNodesByX(std::vector<Node>& nodes)
{
    std::sort(nodes.begin(), nodes.end(), 
              [](const Node& a, const Node& b) {
                  return a.x < b.x;
              });
}
#endif