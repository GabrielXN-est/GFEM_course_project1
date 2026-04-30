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

// Não chamar FEM p-hieraquico como PoU do GFEM, pois não é PoU e não foram implementado // Pode chama-lo se não se foi implementado nenhuma função de enriquecimento

// dof order in the nodes (PoU dofs -> node aplyed enrichments ->
    // -> pGFEM generated polinomial enrichments -> pGFEM generated non-polinomial enrichments)

void plot_series(const plotting_data& data, const std::string& title = "", const std::string& path = "./")
{
    matplot::plot(data.x_values, data.u_values);
    matplot::xlabel("x");
    matplot::ylabel("u(x)");;
    matplot::save(path + title + ".png");
}

template <typename T, typename G>
void plot_error(std::vector<T>& x1, std::vector<T>& x2, std::vector<T>& x3,
     std::vector<G>& y1, std::vector<G>& y2, std::vector<G>& y3,
     std::string_view x_label, std::string title,
     std::string path = "./", bool cond_plot = false)
{
    //matplot::loglog(nelem_L, error_FEM, "-o");
    auto f = matplot::figure(true);
    matplot::loglog(x1, y1,"-o")->display_name("FEM");
    matplot::hold(matplot::on);
    matplot::loglog(x2, y2,"-o")->display_name("GFEM");
    matplot::hold(matplot::on);
    matplot::loglog(x3, y3,"-o")->display_name("S-GFEM");
    matplot::hold(matplot::off);
    auto lgd = matplot::legend();
    lgd->location(matplot::legend::general_alignment::bottomleft);
    matplot::xlabel(x_label);
    if (cond_plot)
        matplot::ylabel("Conditioning number");
    else
        matplot::ylabel("Relative error in energy norm");
    matplot::title(title);
    matplot::save(path + title + ".png");
}

void simulation(std::string path, const std::vector<int>& n_elem_L, // n de elementos a iterar
    std::vector<double>& error_reference, std::vector<double>& dofs_reference, std::vector<double>& conditioning_reference, // outputs
    double L, double x_gamma,std::vector<double> E, double A, int bf_func, double U_exact, // parametros fixos
    std::string title, int porder, std::string eltype, std::vector<double> geom_enr = {}, double tol = std::pow(10, -30))
{
    // rerve memory for outputs
    error_reference.reserve(static_cast<int>(n_elem_L.size()));
    dofs_reference.reserve(static_cast<int>(n_elem_L.size()));
    conditioning_reference.reserve(static_cast<int>(n_elem_L.size()));

    for (int nelem: n_elem_L)
    {
        std::string filename {path + "/input_files/" + title + "_" + std::to_string(nelem) + ".txt"};

        // create input file
        if (eltype == "lBar")
            generate_input(filename, nelem, porder, eltype, L, E, {}, A, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
            std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
            std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
            bf_func, 0, 0, 0, x_gamma); // bf_func_id, alpha, xb, xi, xgamma
        else
            generate_input(filename, nelem, 1, eltype, L, E, {}, A, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
            std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
            std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
            bf_func, 0, 0, 0, x_gamma, porder-1, geom_enr); // bf_func_id, alpha, xb, xi, xgamma
        
        // start mesh and compute local constants
        Mesh mesh {};
        read_input(filename, mesh);

        // Assemble global system
        mesh.assemble_direct();

        // Scale the global matrixes and vectors to improve conditioning
        mesh.create_scaled_global_system(true);

        //Solve the system
        if (eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M" || eltype == "pGFEMBar" || eltype == "pGFEMBar_sc" || eltype == "pGFEMBar_2Proj" || eltype == "pSGFEMBar_2Proj")
        {
            try
            {
                mesh.solve_dependent_system(tol, 100000000);
            }
            catch(const std::exception& e)
            {
                std::cout << '\t' << e.what() << std::endl;
                std::cout << "\t Aumentando a tolerancia para " << std::pow(10, -20) << std::endl;
                try
                {mesh.solve_dependent_system(std::pow(10, -20), 1000000);}
                catch(const std::exception& e)
                {
                    std::cout << '\t' << e.what() << std::endl;
                    std::cout << "\t Aumentando a tolerancia para " << std::pow(10, -12) << std::endl;
                    mesh.solve_dependent_system(std::pow(10, -12), 1000000);
                }
            }
        }
        else
            mesh.solve();

        // post processing
        mesh.complete_U();
        error_reference.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
        conditioning_reference.push_back(mesh.scaled_condition_number);
        dofs_reference.push_back(mesh.K_global_pos.mat.size());

        std::cout << "Relative error in energy norm for " << title << " with " << nelem << " elements equals: " << error_reference.back() << std::endl;
        std::cout << "Conditioning number for " << title << " with " << nelem << " elements equals: " << conditioning_reference.back() << std::endl << std::endl;

        plot_series(get_solution_plotable(mesh, L/100), title + "_" + std::to_string(nelem), path + "/plots/");
    }
}

#define P_VERSION
//#define H_VERSION
int main()
{
    try {
        std::string path {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 3/temp/"};
        
        // parameters for the problem
        double L {1};
        double x_gamma {0};
        std::vector<double> E {1.};
        std::vector<double> E_xlim {};
        double A {1};
        int bf_func{20};
        double U_exact {6.425951283957258};

        // solution vectors
        std::vector<double> h_FEM_error {};
        std::vector<double> h_GFEM_topo_error {}, h_GFEM_geom_error {};
        std::vector<double> h_SGFEM_topo_error {}, h_SGFEM_geom_error {};

        // dofs vectors
        std::vector<double> h_FEM_dofs {};
        std::vector<double> h_GFEM_topo_dofs {}, h_GFEM_geom_dofs {};
        std::vector<double> h_SGFEM_topo_dofs {}, h_SGFEM_geom_dofs {};

        // conditioning vectors
        std::vector<double> h_FEM_conditioning {};
        std::vector<double> h_GFEM_topo_conditioning {}, h_GFEM_geom_conditioning {};
        std::vector<double> h_SGFEM_topo_conditioning {}, h_SGFEM_geom_conditioning {};
        
        // number of elements
        std::vector<int> nelem_L{8,16, 32, 64, 128};
    
        std::cout << "________________h-version SGFEM geometrical________________" << std::endl;
            simulation(path, nelem_L, h_SGFEM_geom_error, h_SGFEM_geom_dofs, h_SGFEM_geom_conditioning,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "P2_SGFEM_geom_pord_" + std::to_string(1), 1, "pSGFEMBar_2Proj", {0,1./4.});

        std::cout << "________________Convergence Rates________________" << std::endl;
            // taxas de convergência
            int size {static_cast<int>(nelem_L.size())};

            std::cout << "Convergence rate for geometrical SGFEM in terms of h: " << (std::log(h_SGFEM_geom_error[size-1])-std::log(h_SGFEM_geom_error[size-2]))/(std::log(L/nelem_L[size-1])-std::log(L/nelem_L[size-2])) << "\n\n";

            std::cout << "\nIn respect to dofs:\n";

            std::cout << "Convergence rate for geometrical SGFEM in terms of dofs: " << -(std::log(h_SGFEM_geom_error[size-1])-std::log(h_SGFEM_geom_error[size-2]))/(std::log(h_SGFEM_geom_dofs[size-1])-std::log(h_SGFEM_geom_dofs[size-2])) << "\n\n";
            std::cout << "\n";

        std::cout << "________________Condition number growth rate________________" << std::endl;
                    std::cout << "\nIn respect to dofs:\n";
            std::cout << "Convergence rate for geometrical SGFEM in terms of dofs: " << (std::log(h_SGFEM_geom_conditioning[size-1])-std::log(h_SGFEM_geom_conditioning[size-2]))/(std::log(h_SGFEM_geom_dofs[size-1])-std::log(h_SGFEM_geom_dofs[size-2])) << "\n\n";
            
        std::cout << "________________Plotagem________________" << std::endl;
            //plotagens
            plot_error(h_FEM_dofs, h_GFEM_geom_dofs, h_SGFEM_geom_dofs, h_FEM_error, h_GFEM_geom_error, h_SGFEM_geom_error, "Degrees of freedom", "Geometrical enrichment error", path + "/plots/convergence/");
            plot_error(h_FEM_dofs, h_GFEM_geom_dofs, h_SGFEM_geom_dofs, h_FEM_conditioning, h_GFEM_geom_conditioning, h_SGFEM_geom_conditioning, "Degrees of freedom", "Geometrical enrichment conditioning", path + "/plots/convergence/", true);
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