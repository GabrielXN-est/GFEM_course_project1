#include "Enrichment.h"

Node::~Node()
{
    for (Enrichment* e : enr)
        delete e;
}

void Node::get_polinomial_order_of_enrichment()
{
    for (Enrichment* e : enr)
    {
        if (e->grau > polinomial_order_of_enrichment)
            polinomial_order_of_enrichment = e->grau;
    }
}

void sortNodesByX(std::vector<Node*>& nodes)
{
    std::sort(nodes.begin(), nodes.end(), 
              [](const Node* a, const Node* b) {
                  return a->x < b->x;
              });
}

void sortNodesByX(std::vector<Node>& nodes)
{
    std::sort(nodes.begin(), nodes.end(), 
              [](const Node& a, const Node& b) {
                  return a.x < b.x;
              });
}