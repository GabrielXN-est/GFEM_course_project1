#include "plot_solution.h"

void plot_series(std::vector<plotting_data>& data, const std::string& title)
{
    for (const plotting_data& datai : data)
        {matplot::plot(datai.x_values, datai.u_values, datai.label);}
    matplot::xlabel("x");
    matplot::ylabel("u(x)");
    matplot::legend();
}

plotting_data get_solution_plotable(Mesh& mesh, double precision = 0.01, std::string label = "")
{
    sortNodesByX(mesh.nodes);
    
    double x_min {mesh.nodes[0].x};
    double x_max {mesh.nodes[mesh.nodes.size()-1].x};

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
                for (std::size_t i {0}; i < el->Nod_list.size(); i++)
                    u_x += N.vec[i]*mesh.U[el->conectivity[i]];
                break;
            }
        }
    }

    return {x_values, u_values, label};
}