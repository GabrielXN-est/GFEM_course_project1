#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <type_traits>
#include <fstream>

#include "create_input.h"
#include "mesh.h"

#include <matplot/matplot.h>

// Não chamar FEM p-hieraquico como PoU do GFEM, pois não é PoU e não foram implementado // Pode chama-lo se não se foi implementado nenhuma função de enriquecimento

// dof order in the nodes (PoU dofs -> node aplyed enrichments ->
    // -> pGFEM generated polinomial enrichments -> pGFEM generated non-polinomial enrichments)
void plot_error(std::vector<int>& nelem_L, std::vector<double>& error_FEM, std::vector<double>& error_GFEM)
{
    matplot::loglog(nelem_L, error_FEM, "-o");
    matplot::loglog(nelem_L, error_GFEM, "-o");
    matplot::xlabel("Number of elements");
    matplot::ylabel("Relative error in energy norm");
    matplot::show();
}

double convergence_rate(std::vector<int>& nelem_L, std::vector<double>& error, double L)
{
    std::vector<double> log_h {}, log_error {};
    for (std::size_t i {0}; i < nelem_L.size(); i++)
    {
        log_h.push_back(std::log(L/nelem_L[i]));
        log_error.push_back(std::log(error[i]));
    }
    double x_ {0}, y_ {0}, Sxy {0}, Sx2 {0};
    for (std::size_t i {0}; i < log_h.size(); i++)
    {   x_ += log_h[i]; 
        y_ += log_error[i];}

    x_ /= log_h.size();
    y_ /= log_error.size();
    
    for (std::size_t i {0}; i < log_h.size(); i++)
    {   Sxy += (log_h[i]-x_)*(log_error[i]-y_);
        Sx2 += (log_h[i]-x_)*(log_h[i]-x_);}

    return (Sxy/Sx2);

}

int main()
{
    try {
        double U_exact {0.0408777548};
        std::vector<int> nelem_L{2, 4, 8, 16, 32};

        std::cout << "________________h-version FEM________________" << std::endl;
        std::vector<double> h_FEM_error;
        h_FEM_error.reserve(5);
        
        /*
        for (int nelem: nelem_L)
        {
            std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_3_FEM_nel_" + std::to_string(nelem) + ".txt"};

            generate_input(filename, nelem, 2, "lBar", 1, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
            std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
            std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
            10, 0.5, 0.2); // bf_func_id, alpha, xb

            Mesh mesh {};
            read_input(filename, mesh);

            mesh.assemble_direct();
            mesh.solve();

            mesh.complete_U();
            h_FEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
            std::cout << "Relative error in energy norm for" << nelem << " elements: " << h_FEM_error.back() << std::endl;
        }
        */
        std::cout << "________________h-version GFEM________________" << std::endl;
        std::vector<double> h_GFEM_error;
        h_GFEM_error.reserve(5);
        for (int nelem: nelem_L)
        {
            std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_3_GFEM_nel_" + std::to_string(nelem) + ".txt"};

            generate_input(filename, nelem, 1, "pGFEMBar", 1, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
            std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
            std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
            10, 0.5, 0.2, 0, 0, 1); // bf_func_id, alpha, xb, xi, xgamma, porder_Enrichment

            Mesh mesh {};
            read_input(filename, mesh);

            mesh.assemble_direct();
            mesh.solve();

            mesh.complete_U();
            h_GFEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
            std::cout << "Relative error in energy norm for" << nelem << " elements: " << h_FEM_error.back() << std::endl;
        }


        plot_error(nelem_L, h_FEM_error, h_GFEM_error);
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Error: " << "undefined" << std::endl;
        return 1;
    }
}