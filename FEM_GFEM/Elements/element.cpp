#include "element.h"

//Element
void Element::get_nodes(std::vector<Node>& nodevec)
{
    int j {0};
    for (Node& i: nodevec)
    {
        if (std::count(N_list.begin(), N_list.end(), i.id))
        {
            i.vicinal_elements.push_back(this);
            Nod_list[j++] = &i;
        }
    }
    sortNodesByX(Nod_list);
}

void Element::Set_ndof()
{
    Ndof = 0;
    for (Node* node : Nod_list)
        {Ndof += static_cast<int>(node->dofs.size());}
    K_local.inicialize(static_cast<std::size_t>(Ndof), static_cast<std::size_t>(Ndof));
    F_local.inicialize(static_cast<std::size_t>(Ndof));
}