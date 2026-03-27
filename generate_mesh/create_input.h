#include <vector>
#include <string>

void generate_input(std::string filename, int nel, int porder, std::string eltype,
    double L, const std::vector<double>& E, const std::vector<double>& Exlim, double A, double C,
    const std::vector<double>& d_bcs, const std::vector<int>& d_bcs_pos, const std::vector<int>& d_bcs_dofs, // dirichilet boundary conditions
    const std::vector<double>& f_bcs, const std::vector<int>& f_bcs_pos, const std::vector<int>& f_bcs_dofs, // Neumann Boundary conditions
    int bf_func_id, double alpha=0.0, double xb=0.0, double xi = 0.0, 
    double xgamma = 0.0, int porder_Enrichment = 1);