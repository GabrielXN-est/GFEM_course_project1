#include "plot_solution.h"

plotting_data get_solution_plotable(Mesh& mesh, double precision, std::string label)
{
    double x_min {mesh.Big_number};
    double x_max {-mesh.Big_number};
    for (Node& node: mesh.nodes)
        if (node.x < x_min)
            x_min = node.x;
        else if (node.x > x_max)
            x_max = node.x;

    std::vector<double> x_values {}, u_values {};  
    
    x_values.reserve(static_cast<std::size_t>((x_max-x_min)/precision)+1);
    u_values.reserve(x_values.size());

    for (double x = x_min; x <= x_max; x += precision)
    {
        x_values.push_back(x);
        double u_x {0};

        for (Element* el : mesh.c_bars)
        {
            if (el->Nod_list[0]->x <= x && el->Nod_list[0]->x + el->el_size >= x)
            {
                shape_functions* sf {el->get_shape_func()};
                sf->operator()(el->mapping(x, el->Nod_list[0]->x, el->el_size));
                sf->mont_vector();
                Vector N (el->Ndof);
                el->Mont_N(N, sf, el->Nod_list, x, el->Ndof);
                for (std::size_t i {0}; i < N.vec.size(); i++)
                    u_x += N.vec[i]*mesh.U[el->conectivity[i]];
                break;
            }
        }
        u_values.push_back(u_x);
    }

    return {x_values, u_values, label};
}