#include "mesh.h"

struct plotting_data
{
    std::vector<double> x_values {};
    std::vector<double> u_values {};
    std::string label {};
};

plotting_data get_solution_plotable(Mesh& mesh, double precision = 0.01, std::string label = "");