#include "mesh.h"
#include <matplot/matplot.h>

struct plotting_data
{
    std::vector<double> x_values {};
    std::vector<double> u_values {};
    std::string label {};
};

void plot_series(std::vector<plotting_data>& data, const std::string& title);

plotting_data get_solution_plotable(Mesh& mesh, double precision = 0.01, std::string label = "");