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
    matplot::plot(data.x_values, data.u_values, data.label);
    matplot::xlabel("x");
    matplot::ylabel("u(x)");
    matplot::legend();
    matplot::save(path + title + ".png");
}
template <typename T, typename G>
void plot_error(std::vector<T>& x1, std::vector<T>& x2, std::vector<G>& y1, std::vector<G>& y2, std::string_view x_label, std::string_view title)
{
    //matplot::loglog(nelem_L, error_FEM, "-o");
    auto f = matplot::figure(true);
    matplot::loglog(x1, y1,"-o")->display_name("h-FEM");
    matplot::hold(matplot::on);
    matplot::loglog(x2, y2,"-o")->display_name("h-GFEM");
    matplot::legend();
    matplot::xlabel(x_label);
    matplot::ylabel("Relative error in energy norm");
    matplot::title(title);
    matplot::show();
}
template <typename T, typename G>
void plot_error(std::vector<T>& x1, std::vector<T>& x2, std::vector<T>& x3, std::vector<T>& x4, 
    std::vector<G>& y1, std::vector<G>& y2, std::vector<G>& y3, std::vector<G>& y4, std::string_view x_label, std::string_view title)
{
    //matplot::loglog(nelem_L, error_FEM, "-o");
    auto f = matplot::figure(true);
    matplot::loglog(x1, y1,"-o")->display_name("h-FEM");
    matplot::hold(matplot::on);
    matplot::loglog(x2, y2,"-o")->display_name("h-GFEM");
    matplot::hold(matplot::on);
    matplot::loglog(x3, y3,"-o")->display_name("p-FEM");
    matplot::hold(matplot::on);
    matplot::loglog(x4, y4,"-o")->display_name("p-GFEM");
    matplot::legend();
    matplot::xlabel(x_label);
    matplot::ylabel("Relative error in energy norm");
    matplot::title(title);
    matplot::show();
}


int main()
{
    try {
        std::vector<int> pord_L {1,2,3,4,5,6};

        double L {1};

        std::cout << "________________alpha = 0,5________________" << std::endl;
        {
            std::vector<int> nelem_L{2, 4, 8, 16, 32};

            double alpha {0.5};
            double U_exact {0.0408777548};
            std::cout << "________________h-version FEM________________" << std::endl;
            std::vector<double> h_FEM_error {}, h_FEM_dofs {};
            h_FEM_error.reserve(5);
            h_FEM_dofs.reserve(5);
            for (int nelem: nelem_L)
            {
                std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_3_FEM_nel_" + std::to_string(nelem) + ".txt"};

                generate_input(filename, nelem, 2, "lBar", L, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
                std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
                std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
                10, alpha, 0.2); // bf_func_id, alpha, xb

                Mesh mesh {};
                read_input(filename, mesh);

                mesh.assemble_direct();
                mesh.solve();

                mesh.complete_U();
                h_FEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                h_FEM_dofs.push_back(mesh.K_global_pos.mat.size());
                std::cout << "Relative error in energy norm for " << nelem << " elements: " << h_FEM_error.back() << std::endl;
                plot_series(get_solution_plotable(mesh), "EX1_3_FEM_nel_" + std::to_string(nelem) + ".txt", "/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/plots/");
            }

            std::cout << "________________h-version GFEM________________" << std::endl;
            std::vector<double> h_GFEM_error {}, h_GFEM_dofs {};
            h_GFEM_error.reserve(5);
            h_GFEM_dofs.reserve(5);
            for (int nelem: nelem_L)
            {
                std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_3_GFEM_nel_" + std::to_string(nelem) + ".txt"};

                generate_input(filename, nelem, 1, "pGFEMBar", 1, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
                std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
                std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
                10, alpha, 0.2, 0, 0, 1); // bf_func_id, alpha, xb, xi, xgamma, porder_Enrichment

                Mesh mesh {};
                read_input(filename, mesh);

                mesh.assemble_direct();
                mesh.solve_dependent_system(std::pow(10, -30), 10000);

                mesh.complete_U();
                h_GFEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                h_GFEM_dofs.push_back(mesh.K_global_pos.mat.size());
                std::cout << "Relative error in energy norm for " << nelem << " elements: " << h_GFEM_error.back() << std::endl;
            }

            std::cout << "________________p-version FEM________________" << std::endl;
            std::vector<double> p_FEM_error {}, p_FEM_dofs {};
            p_FEM_error.reserve(6);
            p_FEM_dofs.reserve(6);
            for (int pord: pord_L)
            {
                std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_4_FEM_pord_" + std::to_string(pord) + ".txt"};

                generate_input(filename, 2, pord, "lBar", L, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
                std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
                std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
                10, alpha, 0.2); // bf_func_id, alpha, xb

                Mesh mesh {};
                read_input(filename, mesh);

                mesh.assemble_direct();
                mesh.solve();

                mesh.complete_U();
                p_FEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                p_FEM_dofs.push_back(mesh.K_global_pos.mat.size());
                std::cout << "Relative error in energy norm for " << pord << " polynomial order: " << p_FEM_error.back() << std::endl;
            }

            std::cout << "________________p-version GFEM________________" << std::endl;
            std::vector<double> p_GFEM_error {}, p_GFEM_dofs {};
            p_GFEM_error.reserve(6);
            p_GFEM_dofs.reserve(6);
            for (int pord: pord_L)
            {
                std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_4_GFEM_pord_" + std::to_string(pord) + ".txt"};

                generate_input(filename, 2, 1, "pGFEMBar", 1, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
                std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
                std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
                10, alpha, 0.2, 0, 0, pord-1); // bf_func_id, alpha, xb, xi, xgamma, porder_Enrichment

                Mesh mesh {};
                read_input(filename, mesh);

                mesh.assemble_direct();
                mesh.solve_dependent_system(std::pow(10, -30), 100000);

                mesh.complete_U();
                p_GFEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                p_GFEM_dofs.push_back(mesh.K_global_pos.mat.size());
                std::cout << "Relative error in energy norm for " << pord << " polynomial order: " << p_GFEM_error.back() << std::endl;
            }

            // taxas de convergência
            int size = nelem_L.size();
            std::cout << std::endl <<"Convergence rate for h-FEM in terms of h: " << (std::log(h_FEM_error[size-1])-std::log(h_FEM_error[size-2]))/(std::log(L/nelem_L[size-1])-std::log(L/nelem_L[size-2])) << "\n";
            std::cout << "Convergence rate for h-GFEM in terms of h: " << (std::log(h_GFEM_error[size-1])-std::log(h_GFEM_error[size-2]))/(std::log(L/nelem_L[size-1])-std::log(L/nelem_L[size-2])) << "\n";
            std::cout << "Convergence rate for h-FEM in terms of dofs: " << (std::log(h_FEM_error[size-1])-std::log(h_FEM_error[size-2]))/(std::log(h_FEM_dofs[size-1])-std::log(h_FEM_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for h-GFEM in terms of dofs: " << (std::log(h_GFEM_error[size-1])-std::log(h_GFEM_error[size-2]))/(std::log(h_GFEM_dofs[size-1])-std::log(h_GFEM_dofs[size-2])) << "\n";
            size = pord_L.size();
            std::cout << "Convergence rate for p-FEM in terms of dofs: " << (std::log(p_FEM_error[size-1])-std::log(p_FEM_error[size-2]))/(std::log(p_FEM_dofs[size-1])-std::log(p_FEM_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for p-GFEM in terms of dofs: " << (std::log(p_GFEM_error[size-1])-std::log(p_GFEM_error[size-2]))/(std::log(p_GFEM_dofs[size-1])-std::log(p_GFEM_dofs[size-2])) << "\n";
            //plotagens
            plot_error(nelem_L, nelem_L, h_FEM_error, h_GFEM_error, "Number of elements", "alpha = 0,5|h-version convergence");
            plot_error(h_FEM_dofs, h_GFEM_dofs, h_FEM_error, h_GFEM_error, "Number of dofs", "alpha = 0,5|h-version convergence");
            plot_error(h_FEM_dofs, h_GFEM_dofs, p_FEM_dofs, p_GFEM_dofs, h_FEM_error, h_GFEM_error, p_FEM_error, p_GFEM_error,"Number of dofs", "alpha = 0,5|p-version convergence");
        }

        std::cout << "________________alpha = 50________________" << std::endl;
        {
            double alpha {50.};
            double U_exact {25.138142063};

            std::vector<int> nelem_L{5, 10, 20, 40};
            std::cout << "________________h-version FEM________________" << std::endl;
            std::vector<double> h_FEM_error {}, h_FEM_dofs {};
            h_FEM_error.reserve(4);
            h_FEM_dofs.reserve(4);
            for (int nelem: nelem_L)
            {
                std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_5_FEM_nel_" + std::to_string(nelem) + ".txt"};

                generate_input(filename, nelem, 2, "lBar", L, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
                std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
                std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
                10, alpha, 0.2); // bf_func_id, alpha, xb

                Mesh mesh {};
                read_input(filename, mesh);

                mesh.assemble_direct();
                mesh.solve();

                mesh.complete_U();
                h_FEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                h_FEM_dofs.push_back(mesh.K_global_pos.mat.size());
                std::cout << "Relative error in energy norm for " << nelem << " elements: " << h_FEM_error.back() << std::endl;
            }

            std::cout << "________________h-version GFEM________________" << std::endl;
            std::vector<double> h_GFEM_error {}, h_GFEM_dofs {};
            h_GFEM_error.reserve(4);
            h_GFEM_dofs.reserve(4);
            for (int nelem: nelem_L)
            {
                std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_5_GFEM_nel_" + std::to_string(nelem) + ".txt"};

                generate_input(filename, nelem, 1, "pGFEMBar", 1, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
                std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
                std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
                10, alpha, 0.2, 0, 0, 1); // bf_func_id, alpha, xb, xi, xgamma, porder_Enrichment

                Mesh mesh {};
                read_input(filename, mesh);

                mesh.assemble_direct();
                mesh.solve_dependent_system(std::pow(10, -30), 10000);

                mesh.complete_U();
                h_GFEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                h_GFEM_dofs.push_back(mesh.K_global_pos.mat.size());
                std::cout << "Relative error in energy norm for " << nelem << " elements: " << h_GFEM_error.back() << std::endl;
            }

            std::cout << "________________p-version FEM________________" << std::endl;
            std::vector<double> p_FEM_error {}, p_FEM_dofs {};
            p_FEM_error.reserve(6);
            p_FEM_dofs.reserve(6);
            for (int pord: pord_L)
            {
                std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_6_FEM_pord_" + std::to_string(pord) + ".txt"};

                generate_input(filename, 2, pord, "lBar", L, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
                std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
                std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
                10, alpha, 0.2); // bf_func_id, alpha, xb

                Mesh mesh {};
                read_input(filename, mesh);

                mesh.assemble_direct();
                mesh.solve();

                mesh.complete_U();
                p_FEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                p_FEM_dofs.push_back(mesh.K_global_pos.mat.size());
                std::cout << "Relative error in energy norm for " << pord << " polynomial order: " << p_FEM_error.back() << std::endl;
            }

            std::cout << "________________p-version GFEM________________" << std::endl;
            std::vector<double> p_GFEM_error {}, p_GFEM_dofs {};
            p_GFEM_error.reserve(6);
            p_GFEM_dofs.reserve(6);
            for (int pord: pord_L)
            {
                std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_files/EX1_6_GFEM_pord_" + std::to_string(pord) + ".txt"};

                generate_input(filename, 2, 1, "pGFEMBar", 1, std::vector<double> {1}, std::vector<double> {}, 1, 0, // filename, nel, porder, eltype, L, E, Exlim, A, C,
                std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<int> {1, 1}, // d_bcs, d_bcs_pos, d_bcs_dofs,
                std::vector<double> {}, std::vector<int> {}, std::vector<int> {},// f_bcs, f_bcs_pos, f_bcs_dofs,
                10, alpha, 0.2, 0, 0, pord-1); // bf_func_id, alpha, xb, xi, xgamma, porder_Enrichment

                Mesh mesh {};
                read_input(filename, mesh);

                mesh.assemble_direct();
                mesh.solve_dependent_system(std::pow(10, -30), 100000);

                mesh.complete_U();
                p_GFEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                p_GFEM_dofs.push_back(mesh.K_global_pos.mat.size());
                std::cout << "Relative error in energy norm for " << pord << " polynomial order: " << p_GFEM_error.back() << std::endl;
            }

            // taxas de convergência
            int size = nelem_L.size();
            std::cout << std::endl <<"Convergence rate for h-FEM in terms of h: " << (std::log(h_FEM_error[size-1])-std::log(h_FEM_error[size-2]))/(std::log(L/nelem_L[size-1])-std::log(L/nelem_L[size-2])) << "\n";
            std::cout << "Convergence rate for h-GFEM in terms of h: " << (std::log(h_GFEM_error[size-1])-std::log(h_GFEM_error[size-2]))/(std::log(L/nelem_L[size-1])-std::log(L/nelem_L[size-2])) << "\n";
            std::cout << "Convergence rate for h-FEM in terms of dofs: " << (std::log(h_FEM_error[size-1])-std::log(h_FEM_error[size-2]))/(std::log(h_FEM_dofs[size-1])-std::log(h_FEM_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for h-GFEM in terms of dofs: " << (std::log(h_GFEM_error[size-1])-std::log(h_GFEM_error[size-2]))/(std::log(h_GFEM_dofs[size-1])-std::log(h_GFEM_dofs[size-2])) << "\n";
            size = pord_L.size();
            std::cout << "Convergence rate for p-FEM in terms of dofs: " << (std::log(p_FEM_error[size-1])-std::log(p_FEM_error[size-2]))/(std::log(p_FEM_dofs[size-1])-std::log(p_FEM_dofs[size-2])) << "\n";
            std::cout << "Convergence rate for p-GFEM in terms of dofs: " << (std::log(p_GFEM_error[size-1])-std::log(p_GFEM_error[size-2]))/(std::log(p_GFEM_dofs[size-1])-std::log(p_GFEM_dofs[size-2])) << "\n";
            //plotagens
            plot_error(nelem_L, nelem_L, h_FEM_error, h_GFEM_error, "Number of elements", "alpha = 50|h-version convergence");
            plot_error(h_FEM_dofs, h_GFEM_dofs, h_FEM_error, h_GFEM_error, "Number of dofs", "alpha = 50|h-version convergence");
            plot_error(h_FEM_dofs, h_GFEM_dofs, p_FEM_dofs, p_GFEM_dofs, h_FEM_error, h_GFEM_error, p_FEM_error, p_GFEM_error,"Number of dofs", "alpha = 50|p-version convergence");
        }
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