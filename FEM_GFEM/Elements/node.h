#ifndef NODE_H
#define NODE_H

#include <vector>
#include <algorithm>

// agrupamento de dofs
class Node
{
    public:
    int id;
    double x {};
    std::vector<int> dofs {};

    Node (int index, double x_coord) : id{index}, x {x_coord}{}
    Node ();
    
};

// Função para ordenar um vetor de nodes pela posição x
inline void sortNodesByX(std::vector<Node*>& nodes)
{
    std::sort(nodes.begin(), nodes.end(), 
              [](const Node* a, const Node* b) {
                  return a->x < b->x;
              });
}

#endif