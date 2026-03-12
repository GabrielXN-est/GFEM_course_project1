#include "element.h"

//Element
void Element::get_nodes(std::vector<Node>& nodevec)
{
    int j {0};
    for (Node& i: nodevec)
    {
        if (std::count(N_list.begin(), N_list.end(), i.id))
            Nod_list[j++] = &i;
    }
    sortNodesByX(Nod_list);
}