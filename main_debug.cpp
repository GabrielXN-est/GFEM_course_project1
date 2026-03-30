#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <type_traits>
#include <fstream>
#include <thread>

#include "create_input.h"
#include "mesh.h"
#include "plot_solution.h"

#include <matplot/matplot.h>

int main()
{
    generate_input("debug.txt", 3, 1, "pGFEMBar_WD_M", 12, {1}, {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
    std::vector<double> {{0.,1.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
    std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
    12, 0, 0, 0, 5); // bf_func_id, alpha, xb, xi, xgamma

    Mesh mesh {};
    read_input("debug.txt", mesh);

    Moes_enrichment_1D enr {1, 5, &mesh.nodes[1]};
    std::vector<double> x_values {}, u_values {}, uD_values {}; 

    for (double x {0}; x <= 8; x += 0.1)
    
    {
        x_values.push_back(x);
        u_values.push_back(enr(x));
        uD_values.push_back(enr.D(x));
    }

    matplot::plot(x_values, u_values);
    matplot::hold(matplot::on);
    matplot::plot(x_values, uD_values);
    matplot::xlabel("x");
    matplot::ylabel("u(x)");
    matplot::show();
    return 0;
}