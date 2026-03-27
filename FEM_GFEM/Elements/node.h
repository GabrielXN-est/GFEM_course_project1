#ifndef NODE_H
#define NODE_H

#include <vector>
#include <algorithm>

class Element;
class Enrichment;

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
    int ep_enriched {-1}; // flag para indicar se o nó já é enriquecido outros enriquecimentos com polinomios
    int polinomial_order_of_enrichment {0}; // indica o grau do polinômio de enriquecimento do nó
    double biggest_vicinal_element_size {}; // para enriquecimentos polinomiais, indica o tamanho do maior elemento vicinal ao nó

    Node (int index, double x_coord, std::vector<int> enrichment_ids) : id{index}, x {x_coord}, enr_ids{enrichment_ids}{}
    Node ();
    
    ~Node ();

    void get_polinomial_order_of_enrichment();
};

// Função para ordenar um vetor de nodes pela posição x
void sortNodesByX(std::vector<Node*>& nodes);

#endif