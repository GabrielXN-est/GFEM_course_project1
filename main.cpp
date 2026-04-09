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
void plot_error(std::vector<T>& x1, std::vector<T>& x2, std::vector<G>& y1, std::vector<G>& y2, std::string_view x_label, std::string title,
     std::string path = "./")
{
    //matplot::loglog(nelem_L, error_FEM, "-o");
    auto f = matplot::figure(true);
    matplot::loglog(x1, y1,"-o")->display_name("pair number of elements");
    matplot::hold(matplot::on);
    matplot::loglog(x2, y2,"-o")->display_name("odd number of elements");
    matplot::hold(matplot::off);
    auto lgd = matplot::legend();
    lgd->location(matplot::legend::general_alignment::bottomleft);
    matplot::xlabel(x_label);
    matplot::ylabel("Relative error in energy norm");
    matplot::title(title);
    matplot::save(path + title + ".png");
}

void simulation(const std::vector<int>& n_elem_L, std::vector<double>& error_reference, std::vector<double>& dofs_reference,
double L, double x_gamma,std::vector<double> E, double A, int bf_func, double U_exact, // parametros fixos
std::string title, int porder, std::string eltype, double tol = std::pow(10, -30))
{
    std::vector<double> E_xlim {x_gamma, L};

    error_reference.reserve(static_cast<int>(n_elem_L.size()));
    dofs_reference.reserve(static_cast<int>(n_elem_L.size()));

    for (int nelem: n_elem_L)
    {
        std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/" + title + "_" + std::to_string(nelem) + ".txt"};

        if (eltype == "lBar")
            generate_input(filename, nelem, porder, eltype, L, E, E_xlim, A, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
            std::vector<double> {{0.,1.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
            std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
            bf_func, 0, 0, 0, x_gamma); // bf_func_id, alpha, xb, xi, xgamma
        else
            generate_input(filename, nelem, 1, eltype, L, E, E_xlim, A, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
            std::vector<double> {{0.,1.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
            std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
            bf_func, 0, 0, 0, x_gamma, porder-1); // bf_func_id, alpha, xb, xi, xgamma
        Mesh mesh {};
        read_input(filename, mesh);

        mesh.assemble_direct();
        if (eltype == "pGFEMBar_WD_S" || eltype == "pGFEMBar_WD_M" || eltype == "pGFEMBar" || eltype == "pGFEMBar_sc")
        {
            try
            {
                mesh.solve_dependent_system(tol, 1000000);
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

        mesh.complete_U();
        error_reference.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
        dofs_reference.push_back(mesh.K_global_pos.mat.size());
        std::cout << "Relative error in energy norm for " << title << " with " << nelem << " elements equals: " << error_reference.back() << std::endl;
        plot_series(get_solution_plotable(mesh, L/100), title + "_" + std::to_string(nelem), "/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/plots/");
    }
}

//#define P_VERSION
#define H_VERSION
int main()
{
    try {
        // parameters for the problem
        double L {10};
        double x_gamma {L/2};
        std::vector<double> E {10000, 1000};
        std::vector<double> E_xlim {x_gamma, L};
        double A {1};
        int bf_func{12};
        double U_exact {687125./7392.};

        // solution vectors
        std::vector<double> h_FEM_lin_p_error {}, h_FEM_lin_o_error {}, h_FEM_quad_p_error {}, h_FEM_quad_o_error {};
        std::vector<double> h_GFEM_lin_p_error {}, h_GFEM_lin_o_error {}, h_GFEM_quad_p_error {}, h_GFEM_quad_o_error {};
        std::vector<double> h_GFEM_lin_p_error_suk {}, h_GFEM_lin_o_error_suk {}, h_GFEM_quad_p_error_suk {}, h_GFEM_quad_o_error_suk {};
        std::vector<double> h_GFEM_lin_p_error_moes {}, h_GFEM_lin_o_error_moes {}, h_GFEM_quad_p_error_moes {}, h_GFEM_quad_o_error_moes {};

        std::vector<double> p_FEM_p_error {}, p_FEM_o_error {};
        std::vector<double> p_GFEM_p_error {}, p_GFEM_o_error {};
        std::vector<double> p_GFEM_p_error_suk {}, p_GFEM_o_error_suk {};
        std::vector<double> p_GFEM_p_error_moes {}, p_GFEM_o_error_moes {};

        // dofs vectors
        std::vector<double> h_FEM_lin_p_dofs {}, h_FEM_lin_o_dofs{}, h_FEM_quad_p_dofs {}, h_FEM_quad_o_dofs {};
        std::vector<double> h_GFEM_lin_p_dofs {}, h_GFEM_lin_o_dofs {}, h_GFEM_quad_p_dofs {}, h_GFEM_quad_o_dofs {};
        std::vector<double> h_GFEM_lin_p_dofs_suk {}, h_GFEM_lin_o_dofs_suk {}, h_GFEM_quad_p_dofs_suk {}, h_GFEM_quad_o_dofs_suk {};
        std::vector<double> h_GFEM_lin_p_dofs_moes {}, h_GFEM_lin_o_dofs_moes {}, h_GFEM_quad_p_dofs_moes {}, h_GFEM_quad_o_dofs_moes {};

        std::vector<double> p_FEM_p_dofs {}, p_FEM_o_dofs{};
        std::vector<double> p_GFEM_p_dofs {}, p_GFEM_o_dofs {};
        std::vector<double> p_GFEM_p_dofs_suk {}, p_GFEM_o_dofs_suk {};
        std::vector<double> p_GFEM_p_dofs_moes {}, p_GFEM_o_dofs_moes {};
        #ifdef H_VERSION
        std::cout << "________________h-version FEM linear________________" << std::endl;

        std::vector<int> nelem_p_L{2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
        std::vector<int> nelem_o_L{3, 5, 7, 9, 11, 13, 15, 17, 19};

            std::cout << " \n -> Pairs" << std::endl;
            simulation(nelem_p_L, h_FEM_lin_p_error, h_FEM_lin_p_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_FEM_pord_" + std::to_string(1), 1, "lBar");

            std::cout << " \n -> Odds" << std::endl;
            simulation(nelem_o_L, h_FEM_lin_o_error, h_FEM_lin_o_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_FEM_pord_" + std::to_string(1), 1, "lBar");
 
        std::cout << "________________h-version FEM quadratic________________" << std::endl;
            std::cout << " \n -> Pairs" << std::endl;
            simulation(nelem_p_L, h_FEM_quad_p_error, h_FEM_quad_p_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_FEM_pord_" + std::to_string(2), 2, "lBar");

            std::cout << " \n -> Odds" << std::endl;
            simulation(nelem_o_L, h_FEM_quad_o_error, h_FEM_quad_o_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_FEM_pord_" + std::to_string(2), 2, "lBar");

        std::cout << "________________h-version GFEM linear________________" << std::endl;
            std::cout << " \n -> Pairs" << std::endl;
            simulation(nelem_p_L, h_GFEM_lin_p_error, h_GFEM_lin_p_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_pord_" + std::to_string(1), 1, "pGFEMBar");

            std::cout << " \n -> Odds" << std::endl;
            simulation(nelem_o_L, h_GFEM_lin_o_error, h_GFEM_lin_o_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_pord_" + std::to_string(1), 1, "pGFEMBar");

        std::cout << "________________h-version GFEM quadratic________________" << std::endl;
            std::cout << " \n -> Pairs" << std::endl;
            simulation(nelem_p_L, h_GFEM_quad_p_error, h_GFEM_quad_p_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_pord_" + std::to_string(2), 2, "pGFEMBar");

            std::cout << " \n -> Odds" << std::endl;
            simulation(nelem_o_L, h_GFEM_quad_o_error, h_GFEM_quad_o_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_pord_" + std::to_string(2), 2, "pGFEMBar");

        std::cout << "________________h-version GFEM linear Sukumar________________" << std::endl;
            std::cout << " \n -> Pairs" << std::endl;
            simulation(nelem_p_L, h_GFEM_lin_p_error_suk, h_GFEM_lin_p_dofs_suk,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_S_pord_" + std::to_string(1), 1, "pGFEMBar_WD_S");

            std::cout << " \n -> Odds" << std::endl;
            simulation(nelem_o_L, h_GFEM_lin_o_error_suk, h_GFEM_lin_o_dofs_suk,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_S_pord_" + std::to_string(1), 1, "pGFEMBar_WD_S");

        std::cout << "________________h-version GFEM quadratic Sukumar________________" << std::endl;
            std::cout << " \n -> Pairs" << std::endl;
            simulation(nelem_p_L, h_GFEM_quad_p_error_suk, h_GFEM_quad_p_dofs_suk,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_S_pord_" + std::to_string(2), 2, "pGFEMBar_WD_S");

            std::cout << " \n -> Odds" << std::endl;
            simulation(nelem_o_L, h_GFEM_quad_o_error_suk, h_GFEM_quad_o_dofs_suk,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_S_pord_" + std::to_string(2), 2, "pGFEMBar_WD_S");

        std::cout << "________________h-version GFEM linear Moes________________" << std::endl;
            std::cout << " \n -> Pairs" << std::endl;
            simulation(nelem_p_L, h_GFEM_lin_p_error_moes, h_GFEM_lin_p_dofs_moes,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_M_pord_" + std::to_string(1), 1, "pGFEMBar_WD_M");

            std::cout << " \n -> Odds" << std::endl;
            simulation(nelem_o_L, h_GFEM_lin_o_error_moes, h_GFEM_lin_o_dofs_moes,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_M_pord_" + std::to_string(1), 1, "pGFEMBar_WD_M");

        std::cout << "________________h-version GFEM quadratic Moes________________" << std::endl;
            std::cout << " \n -> Pairs" << std::endl;
            simulation(nelem_p_L, h_GFEM_quad_p_error_moes, h_GFEM_quad_p_dofs_moes,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_M_pord_" + std::to_string(2), 2, "pGFEMBar_WD_M");

            std::cout << " \n -> Odds" << std::endl;
            simulation(nelem_o_L, h_GFEM_quad_o_error_moes, h_GFEM_quad_o_dofs_moes,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_1_GFEM_M_pord_" + std::to_string(2), 2, "pGFEMBar_WD_M");

        std::cout << "________________Taxas de Convergência (h-version)________________" << std::endl;
            // taxas de convergência
            int size {static_cast<int>(nelem_p_L.size())}, sizeO {static_cast<int>(nelem_o_L.size())};;
            std::cout << "Convergence rate for pair linear h-FEM in terms of h: " << (std::log(h_FEM_lin_p_error[size-1])-std::log(h_FEM_lin_p_error[size-2]))/(std::log(L/nelem_p_L[size-1])-std::log(L/nelem_p_L[size-2])) << "\n";
            std::cout << "Convergence rate for odd linear h-FEM in terms of h: " << (std::log(h_FEM_lin_o_error[sizeO-1])-std::log(h_FEM_lin_o_error[sizeO-2]))/(std::log(L/nelem_o_L[sizeO-1])-std::log(L/nelem_o_L[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair quadratic h-FEM in terms of h: " << (std::log(h_FEM_quad_p_error[size-1])-std::log(h_FEM_quad_p_error[size-2]))/(std::log(L/nelem_p_L[size-1])-std::log(L/nelem_p_L[size-2])) << "\n";
            std::cout << "Convergence rate for odd quadratic h-FEM in terms of h: " << (std::log(h_FEM_quad_o_error[sizeO-1])-std::log(h_FEM_quad_o_error[sizeO-2]))/(std::log(L/nelem_o_L[sizeO-1])-std::log(L/nelem_o_L[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair linear h-GFEM in terms of h: " << (std::log(h_GFEM_lin_p_error[size-1])-std::log(h_GFEM_lin_p_error[size-2]))/(std::log(L/nelem_p_L[size-1])-std::log(L/nelem_p_L[size-2])) << "\n";
            std::cout << "Convergence rate for odd linear h-GFEM in terms of h: " << (std::log(h_GFEM_lin_o_error[sizeO-1])-std::log(h_GFEM_lin_o_error[sizeO-2]))/(std::log(L/nelem_o_L[sizeO-1])-std::log(L/nelem_o_L[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair quadratic h-GFEM in terms of h: " << (std::log(h_GFEM_quad_p_error[size-1])-std::log(h_GFEM_quad_p_error[size-2]))/(std::log(L/nelem_p_L[size-1])-std::log(L/nelem_p_L[size-2])) << "\n";
            std::cout << "Convergence rate for odd quadratic h-GFEM in terms of h: " << (std::log(h_GFEM_quad_o_error[sizeO-1])-std::log(h_GFEM_quad_o_error[sizeO-2]))/(std::log(L/nelem_o_L[sizeO-1])-std::log(L/nelem_o_L[sizeO-2])) << "\n\n";
            
            std::cout << "Convergence rate for pair linear h-GFEM in terms of h with Sukumar: " << (std::log(h_GFEM_lin_p_error_suk[size-1])-std::log(h_GFEM_lin_p_error_suk[size-2]))/(std::log(L/nelem_p_L[size-1])-std::log(L/nelem_p_L[size-2])) << "\n";
            std::cout << "Convergence rate for odd linear h-GFEM in terms of h with Sukumar: " << (std::log(h_GFEM_lin_o_error_suk[sizeO-1])-std::log(h_GFEM_lin_o_error_suk[sizeO-2]))/(std::log(L/nelem_o_L[sizeO-1])-std::log(L/nelem_o_L[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair quadratic h-GFEM in terms of h with Sukumar: " << (std::log(h_GFEM_quad_p_error_suk[size-1])-std::log(h_GFEM_quad_p_error_suk[size-2]))/(std::log(L/nelem_p_L[size-1])-std::log(L/nelem_p_L[size-2])) << "\n";
            std::cout << "Convergence rate for odd quadratic h-GFEM in terms of h with Sukumar: " << (std::log(h_GFEM_quad_o_error_suk[sizeO-1])-std::log(h_GFEM_quad_o_error_suk[sizeO-2]))/(std::log(L/nelem_o_L[sizeO-1])-std::log(L/nelem_o_L[sizeO-2])) << "\n\n";
            
            std::cout << "Convergence rate for pair linear h-GFEM in terms of h with Moes: " << (std::log(h_GFEM_lin_p_error_moes[size-1])-std::log(h_GFEM_lin_p_error_moes[size-2]))/(std::log(L/nelem_p_L[size-1])-std::log(L/nelem_p_L[size-2])) << "\n";
            std::cout << "Convergence rate for odd linear h-GFEM in terms of h with Moes: " << (std::log(h_GFEM_lin_o_error_moes[sizeO-1])-std::log(h_GFEM_lin_o_error_moes[sizeO-2]))/(std::log(L/nelem_o_L[sizeO-1])-std::log(L/nelem_o_L[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair quadratic h-GFEM in terms of h with Moes: " << (std::log(h_GFEM_quad_p_error_moes[size-1])-std::log(h_GFEM_quad_p_error_moes[size-2]))/(std::log(L/nelem_p_L[size-1])-std::log(L/nelem_p_L[size-2])) << "\n";
            std::cout << "Convergence rate for odd quadratic h-GFEM in terms of h with Moes: " << (std::log(h_GFEM_quad_o_error_moes[sizeO-1])-std::log(h_GFEM_quad_o_error_moes[sizeO-2]))/(std::log(L/nelem_o_L[sizeO-1])-std::log(L/nelem_o_L[sizeO-2])) << "\n\n";

            std::cout << "\nEm termos de dofs:\n";

            std::cout << "Convergence rate for pair linear h-FEM in terms of dofs: " << -(std::log(h_FEM_lin_p_error[size-1])-std::log(h_FEM_lin_p_error[size-2]))/(std::log(h_FEM_lin_p_dofs[size-1])-std::log(h_FEM_lin_p_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for odd linear h-FEM in terms of dofs: " << -(std::log(h_FEM_lin_o_error[sizeO-1])-std::log(h_FEM_lin_o_error[sizeO-2]))/(std::log(h_FEM_lin_o_dofs[sizeO-1])-std::log(h_FEM_lin_o_dofs[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair quadratic h-FEM in terms of dofs: " << -(std::log(h_FEM_quad_p_error[size-1])-std::log(h_FEM_quad_p_error[size-2]))/(std::log(h_FEM_quad_p_dofs[size-1])-std::log(h_FEM_quad_p_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for odd quadratic h-FEM in terms of dofs: " << -(std::log(h_FEM_quad_o_error[sizeO-1])-std::log(h_FEM_quad_o_error[sizeO-2]))/(std::log(h_FEM_quad_o_dofs[sizeO-1])-std::log(h_FEM_quad_o_dofs[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair linear h-GFEM in terms of dofs: " << -(std::log(h_GFEM_lin_p_error[size-1])-std::log(h_GFEM_lin_p_error[size-2]))/(std::log(h_GFEM_lin_p_dofs[size-1])-std::log(h_GFEM_lin_p_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for odd linear h-GFEM in terms of dofs: " << -(std::log(h_GFEM_lin_o_error[sizeO-1])-std::log(h_GFEM_lin_o_error[sizeO-2]))/(std::log(h_GFEM_lin_o_dofs[sizeO-1])-std::log(h_GFEM_lin_o_dofs[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair quadratic h-GFEM in terms of dofs: " << -(std::log(h_GFEM_quad_p_error[size-1])-std::log(h_GFEM_quad_p_error[size-2]))/(std::log(h_GFEM_quad_p_dofs[size-1])-std::log(h_GFEM_quad_p_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for odd quadratic h-GFEM in terms of dofs: " << -(std::log(h_GFEM_quad_o_error[sizeO-1])-std::log(h_GFEM_quad_o_error[sizeO-2]))/(std::log(h_GFEM_quad_o_dofs[sizeO-1])-std::log(h_GFEM_quad_o_dofs[sizeO-2])) << "\n\n";
            
            std::cout << "Convergence rate for pair linear h-GFEM in terms of dofs with Sukumar: " << -(std::log(h_GFEM_lin_p_error_suk[size-1])-std::log(h_GFEM_lin_p_error_suk[size-2]))/(std::log(h_GFEM_lin_p_dofs[size-1])-std::log(h_GFEM_lin_p_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for odd linear h-GFEM in terms of dofs with Sukumar: " << -(std::log(h_GFEM_lin_o_error_suk[sizeO-1])-std::log(h_GFEM_lin_o_error_suk[sizeO-2]))/(std::log(h_GFEM_lin_o_dofs[sizeO-1])-std::log(h_GFEM_lin_o_dofs[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair quadratic h-GFEM in terms of dofs with Sukumar: " << -(std::log(h_GFEM_quad_p_error_suk[size-1])-std::log(h_GFEM_quad_p_error_suk[size-2]))/(std::log(h_GFEM_quad_p_dofs[size-1])-std::log(h_GFEM_quad_p_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for odd quadratic h-GFEM in terms of dofs with Sukumar: " << -(std::log(h_GFEM_quad_o_error_suk[sizeO-1])-std::log(h_GFEM_quad_o_error_suk[sizeO-2]))/(std::log(h_GFEM_quad_o_dofs[sizeO-1])-std::log(h_GFEM_quad_o_dofs[sizeO-2])) << "\n\n";
            
            std::cout << "Convergence rate for pair linear h-GFEM in terms of dofs with Moes: " << -(std::log(h_GFEM_lin_p_error_moes[size-1])-std::log(h_GFEM_lin_p_error_moes[size-2]))/(std::log(h_GFEM_lin_p_dofs[size-1])-std::log(h_GFEM_lin_p_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for odd linear h-GFEM in terms of dofs with Moes: " << -(std::log(h_GFEM_lin_o_error_moes[sizeO-1])-std::log(h_GFEM_lin_o_error_moes[sizeO-2]))/(std::log(h_GFEM_lin_o_dofs[sizeO-1])-std::log(h_GFEM_lin_o_dofs[sizeO-2])) << "\n\n";

            std::cout << "Convergence rate for pair quadratic h-GFEM in terms of dofs with Moes: " << -(std::log(h_GFEM_quad_p_error_moes[size-1])-std::log(h_GFEM_quad_p_error_moes[size-2]))/(std::log(h_GFEM_quad_p_dofs[size-1])-std::log(h_GFEM_quad_p_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for odd quadratic h-GFEM in terms of dofs with Moes: " << -(std::log(h_GFEM_quad_o_error_moes[sizeO-1])-std::log(h_GFEM_quad_o_error_moes[sizeO-2]))/(std::log(h_GFEM_quad_o_dofs[sizeO-1])-std::log(h_GFEM_quad_o_dofs[sizeO-2])) << "\n\n";
            
            //plotagens
            plot_error(h_FEM_lin_p_dofs, h_FEM_lin_o_dofs, h_FEM_lin_p_error, h_FEM_lin_o_error, "Degrees of freedom", "Ex2 - Standart linear FEM", "./plots/convergence/");
            plot_error(h_FEM_quad_p_dofs, h_FEM_quad_o_dofs, h_FEM_quad_p_error, h_FEM_quad_o_error, "Degrees of freedom", "Ex2 - Standart quadratic FEM", "./plots/convergence/");
            plot_error(h_GFEM_lin_p_dofs, h_GFEM_lin_o_dofs, h_GFEM_lin_p_error, h_GFEM_lin_o_error, "Degrees of freedom", "Ex2 - linear GFEM = FEM", "./plots/convergence/");
            plot_error(h_GFEM_quad_p_dofs, h_GFEM_quad_o_dofs, h_GFEM_quad_p_error, h_GFEM_quad_o_error, "Degrees of freedom", "Ex2 - quadratic GFEM", "./plots/convergence/");
            plot_error(h_GFEM_lin_p_dofs_suk, h_GFEM_lin_o_dofs_suk, h_GFEM_lin_p_error_suk, h_GFEM_lin_o_error_suk, "Degrees of freedom", "Ex2 - linear GFEM with Sukumar", "./plots/convergence/");
            plot_error(h_GFEM_quad_p_dofs_suk, h_GFEM_quad_o_dofs_suk, h_GFEM_quad_p_error_suk, h_GFEM_quad_o_error_suk, "Degrees of freedom", "Ex2 - quadratic GFEM with Sukumar", "./plots/convergence/");
            plot_error(h_GFEM_lin_p_dofs_moes, h_GFEM_lin_o_dofs_moes, h_GFEM_lin_p_error_moes, h_GFEM_lin_o_error_moes, "Degrees of freedom", "Ex2 - linear GFEM with Moes", "./plots/convergence/");
            plot_error(h_GFEM_quad_p_dofs_moes, h_GFEM_quad_o_dofs_moes, h_GFEM_quad_p_error_moes, h_GFEM_quad_o_error_moes, "Degrees of freedom", "Ex2 - quadratic GFEM with Moes", "./plots/convergence/");
        #endif
        #ifdef P_VERSION
        std::cout << "________________p-version FEM________________" << std::endl;

        std::vector<int> pord_L{1, 2, 3, 4, 5};
        double tol {std::pow(10, -18)};
        std::cout << " \n -> 4 elementos" << std::endl;
        for (int pord: pord_L)
        {
            simulation({4}, p_FEM_p_error, p_FEM_p_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_2_FEM_pord_" + std::to_string(pord), pord, "lBar", tol);
        }
        std::cout << " \n -> 5 elementos" << std::endl;
        for (int pord: pord_L)
        {
            simulation({5}, p_FEM_o_error, p_FEM_o_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_2_FEM_pord_" + std::to_string(pord), pord, "lBar", tol);
        }

        std::cout << "________________p-version GFEM________________" << std::endl;
        std::cout << " \n -> 4 elementos" << std::endl;
        for (int pord: pord_L)
        {
            simulation({4}, p_GFEM_p_error, p_GFEM_p_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_2_GFEM_pord_" + std::to_string(pord), pord, "pGFEMBar", tol);
        }
        std::cout << " \n -> 5 elementos" << std::endl;
        for (int pord: pord_L)
        {
            simulation({5}, p_GFEM_o_error, p_GFEM_o_dofs,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_2_GFEM_pord_" + std::to_string(pord), pord, "pGFEMBar", tol);
        }

        std::cout << "________________p-version GFEM Sukumar________________" << std::endl;
        std::cout << " \n -> 4 elementos" << std::endl;
        for (int pord: pord_L)
        {
            simulation({4}, p_GFEM_p_error_suk, p_GFEM_p_dofs_suk,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_2_GFEM_S_pord_" + std::to_string(pord), pord, "pGFEMBar_WD_S", tol);
        }
        std::cout << " \n -> 5 elementos" << std::endl;
        for (int pord: pord_L)
        {
            simulation({5}, p_GFEM_o_error_suk, p_GFEM_o_dofs_suk,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_2_GFEM_S_pord_" + std::to_string(pord), pord, "pGFEMBar_WD_S", tol);
        }

        std::cout << "________________p-version GFEM Moes________________" << std::endl;
        std::cout << " \n -> 4 elementos" << std::endl;
        for (int pord: pord_L)
        {
            simulation({4}, p_GFEM_p_error_moes, p_GFEM_p_dofs_moes,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_2_GFEM_M_pord_" + std::to_string(pord), pord, "pGFEMBar_WD_M", tol);
        }
        std::cout << " \n -> 5 elementos" << std::endl;
        for (int pord: pord_L)
        {
            simulation({5}, p_GFEM_o_error_moes, p_GFEM_o_dofs_moes,
            L, x_gamma, E, A, bf_func, U_exact, // parametros fixos
            "EX2_2_GFEM_M_pord_" + std::to_string(pord), pord, "pGFEMBar_WD_M", tol);
        }

        //plotagens
        std::cout << "________________Plotagens (p-version)________________" << std::endl;
            plot_error(p_FEM_p_dofs, p_FEM_o_dofs, p_FEM_p_error, p_FEM_o_error, "Degrees of freedom", "Ex2 - Standart FEM (p-version)", "./plots/convergence/");
            plot_error(p_GFEM_p_dofs, p_GFEM_o_dofs, p_GFEM_p_error, p_GFEM_o_error, "Degrees of freedom", "Ex2 - p-GFEM (p-version)", "./plots/convergence/");
            plot_error(p_GFEM_p_dofs_suk, p_GFEM_o_dofs_suk, p_GFEM_p_error_suk, p_GFEM_o_error_suk, "Degrees of freedom", "Ex2 - p-GFEM with Sukumar (p-version)", "./plots/convergence/");
            plot_error(p_GFEM_p_dofs_moes, p_GFEM_o_dofs_moes, p_GFEM_p_error_moes, p_GFEM_o_error_moes, "Degrees of freedom", "Ex2 - p-GFEM with Moes (p-version)", "./plots/convergence/");
        #endif
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