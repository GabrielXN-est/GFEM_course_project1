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

#define FEM
#define GFEM

// Não chamar FEM p-hieraquico como PoU do GFEM, pois não é PoU e não foram implementado // Pode chama-lo se não se foi implementado nenhuma função de enriquecimento

// dof order in the nodes (PoU dofs -> node aplyed enrichments ->
    // -> pGFEM generated polinomial enrichments -> pGFEM generated non-polinomial enrichments)

double u_exact(double x, double Eh, double E1, double E2)
{
    if (x <= 3./8.)
        return x/(Eh);
    else if (x >= 5./8.)
        return 1/Eh*(3./8. + x - 5./8.) + (1./E1+1./E2)/8.;
    else if (static_cast<int>((x-3./8.)*128.) % 2 == 0)
        return 3./8./Eh + 
            1./E1*(x-3./8.-static_cast<double>(static_cast<int>((x-3./8.)*64.))/128.)+
            1./E2*(static_cast<double>(static_cast<int>((x-3./8.)*64.))/128.);
    else
        return 3./8./Eh + 
            1./E2*(x-3./8.-1./128.-static_cast<double>(static_cast<int>((x-3./8.-1./128.)*64.))/128.)+
            1./E1*(1./128.+static_cast<double>(static_cast<int>((x-3./8.-1./128.)*64.))/128.);
}

template <typename T, typename G>
void plot_error(std::vector<T>& x1, std::vector<T>& x2, std::vector<T>& x3, std::vector<T>& x4,
     std::vector<G>& y1, std::vector<G>& y2, std::vector<G>& y3, std::vector<G>& y4,
     std::string_view x_label, std::string title,
     std::string path = "./")
{
    //matplot::loglog(nelem_L, error_FEM, "-o");
    auto f = matplot::figure(true);
    matplot::loglog(x1, y1,"-o")->display_name("FEM");
    matplot::hold(matplot::on);
    matplot::loglog(x2, y2,"-o")->display_name("GFEM-GL - h = L/64");
    matplot::hold(matplot::on);
    matplot::loglog(x3, y3,"-o")->display_name("GFEM-GL - h = L/128");
    matplot::hold(matplot::on);
    matplot::loglog(x4, y4,"-o")->display_name("GFEM-GL - h = L/256");
    matplot::hold(matplot::off);
    auto lgd = matplot::legend();
    lgd->location(matplot::legend::general_alignment::bottomleft);
    matplot::xlabel(x_label);
    matplot::ylabel("Relative error in energy norm");
    matplot::title(title);
    matplot::save(path + title + ".png");
}

int main()
{
    try {
        std::string path {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 3"};
        
        // parameters for the problem
        double L {1.}, A {1.}, P {1.},E1 {1.}, E2 {40.};
        double Eh {E1*E2/(E1*0.5+E2*0.5)};
        int bf_func{0};

        std::vector<double> d_bcs {0.}, f_bcs {P};
        std::vector<int> d_bcs_pos {0}, f_bcs_pos {1};
        std::vector<int> d_bcs_dofs {1}, f_bcs_dofs {1};

        // E
        std::vector<double> E {Eh};
        std::vector<double> E_xlim {3.*L/8.};
        double temp_x {3.*L/8.};
        for (int i {0}; i < 32; i++)
        {
            temp_x += L/128.;
            E_xlim.push_back(temp_x);
            if (i % 2 == 0)
                E.push_back(E1);
            else
                E.push_back(E2);
        }
        E.push_back(Eh);
        E_xlim.push_back(L);

        double U_exact {0.25625};

        // solution vectors
        std::vector<double> h_FEM_error {};
        std::vector<double> h_GFEM_GL_h_L_64_error {}, h_GFEM_GL_h_L_128_error {}, h_GFEM_GL_h_L_256_error {};
        std::vector<std::vector<double>*> h_GFEM_GL_err {&h_GFEM_GL_h_L_64_error, &h_GFEM_GL_h_L_128_error, &h_GFEM_GL_h_L_256_error};

        // dofs vectors
        std::vector<int> h_FEM_dofs{};
        std::vector<int> h_GFEM_GL_h_L_64_dofs {}, h_GFEM_GL_h_L_128_dofs {}, h_GFEM_GL_h_L_256_dofs {};
        std::vector<std::vector<int>*> h_GFEM_GL_dof {&h_GFEM_GL_h_L_64_dofs, &h_GFEM_GL_h_L_128_dofs, &h_GFEM_GL_h_L_256_dofs};

        
        // number of elements
        std::vector<int> nelem_L{8, 16, 32, 64, 128, 256};
        std::vector<int> nelem_L_GL{32, 64, 128};

        // element size
        std::vector<double> H {}, h{};
        H.reserve(nelem_L.size());
        h.reserve(nelem_L_GL.size());

        for (int nelem : nelem_L)
            H.push_back(L/nelem);
        for (int nelem : nelem_L_GL)
            h.push_back(L/nelem);

        // variables to save solutions plotable values for H =1/8
        plotting_data plotable_FEM {};
        std::vector<plotting_data> plotable_GFEM_GL(3);
        # ifdef FEM
        std::cout << "________________h-version FEM________________" << std::endl;
            for (int nelem : nelem_L)
            {
                Mesh mesh {path, "3P_FEM_nelem_" + std::to_string(nelem), nelem,
                    1, "lBar", L, E, E_xlim, A, 0,
                    d_bcs,d_bcs_pos, d_bcs_dofs, // dirichilet boundary conditions
                    f_bcs, f_bcs_pos, f_bcs_dofs, // Neumann Boundary conditions
                    bf_func};
                
                mesh.run();

                h_FEM_error.push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                h_FEM_dofs.push_back(mesh.K_global_pos.mat.size());

                std::cout << "Relative error in energy norm for FEM with " << nelem << " elements equals: " << h_FEM_error.back() << std::endl;

                if (nelem == 8)
                    plotable_FEM = get_solution_plotable(mesh, L/128);
            }
        # endif
        # ifdef GFEM
        int nelem_GL {};
        for (int i = 0; i < nelem_L_GL.size(); i++)
        {
            nelem_GL = nelem_L_GL[i];
            std::cout << "________________h-version GL-GFEM with h = " << L/nelem_GL << "________________" << std::endl;
            for (int nelem : nelem_L)
            {
                Mesh mesh {path, "3P_GFEM_nelem_" + std::to_string(nelem) + "_nelem_GL_" + std::to_string(nelem_GL), nelem,
                    1, "lBar", L, E, E_xlim, A, 0,
                    d_bcs,d_bcs_pos, d_bcs_dofs, // dirichilet boundary conditions
                    f_bcs, f_bcs_pos, f_bcs_dofs, // Neumann Boundary conditions
                    bf_func, 0.0, 0.0, 0.0, 0.0, 0,
                    {3*L/8, 5*L/8}, 1};
                
                mesh.first_run();

                mesh.create_local_problem(L/4, L/2, nelem_GL);
                mesh.run_local_problems();

                mesh.run();

                h_GFEM_GL_err[i]->push_back(std::sqrt(std::abs(U_exact-mesh.strain_energy())/U_exact));
                h_GFEM_GL_dof[i]->push_back(mesh.K_global_pos.mat.size());

                std::cout << "Relative error in energy norm for GFEM with " << nelem << " elements and " << nelem_GL << " local elementos equals: " << h_GFEM_GL_err[i]->back() << std::endl;
            
                if (nelem == 8)
                    plotable_GFEM_GL[i] = get_solution_plotable(mesh, L/128);
            }
        }
        # endif
        std::cout << "________________Convergence Rates________________" << std::endl;
            // taxas de convergência
            int size {static_cast<int>(nelem_L.size())};
            std::cout << "Convergence rate for FEM in terms of dofs: " << -(std::log(h_FEM_error[size-1])-std::log(h_FEM_error[size-2]))/(std::log(h_FEM_dofs[size-1])-std::log(h_FEM_dofs[size-2])) << "\n\n";

            for (int i = 0; i < nelem_L_GL.size(); i++)
            std::cout << "Convergence rate for topological GFEM with h = " << nelem_L_GL[i] << " in terms of dofs: " << 
            (std::log(h_GFEM_GL_err[i]->operator[](size-1))-std::log(h_GFEM_GL_err[i]->operator[](size-2)))/
            (std::log(h_GFEM_GL_dof[i]->operator[](size-1))-std::log(h_GFEM_GL_dof[i]->operator[](size-2))) << "\n";
               
        std::cout << "________________Plotagem________________" << std::endl;
            // get exact solution plotable values
            std::vector<double> x_values {}, u_ex_values {};
            for (double x {0}; x <= L; x += L/1000.)
            {
                x_values.push_back(x);
                u_ex_values.push_back(u_exact(x, Eh, E1, E2));
            }

            //plotagens
            plot_error(h_FEM_dofs, h_GFEM_GL_h_L_64_dofs, h_GFEM_GL_h_L_128_dofs, h_GFEM_GL_h_L_256_dofs,
                       h_FEM_error, h_GFEM_GL_h_L_64_error, h_GFEM_GL_h_L_128_error, h_GFEM_GL_h_L_256_error,
                       "Degrees of freedom", "Error in respect to degrees of freedom", path + "/plots/convergence/");
            plot_error(nelem_L, nelem_L, nelem_L, nelem_L,
                       h_FEM_error, h_GFEM_GL_h_L_64_error, h_GFEM_GL_h_L_128_error, h_GFEM_GL_h_L_256_error,
                       "Global element size", "Error in respect to element size", path + "/plots/convergence/");

            auto f = matplot::figure(true);
            matplot::hold(matplot::on);
            matplot::plot(plotable_FEM.x_values, plotable_FEM.u_values)->display_name("FEM");
            matplot::plot(x_values, u_ex_values)->display_name("Exact solution");
            matplot::plot(plotable_GFEM_GL[0].x_values, plotable_GFEM_GL[0].u_values)->display_name("GFEM GL h=L/64");
            matplot::plot(plotable_GFEM_GL[1].x_values, plotable_GFEM_GL[1].u_values)->display_name("GFEM GL h=L/128");
            matplot::plot(plotable_GFEM_GL[2].x_values, plotable_GFEM_GL[2].u_values)->display_name("GFEM GL h=L/256");
            matplot::hold(matplot::off);
            matplot::xlabel("x");
            matplot::ylabel("u(x)");
            matplot::legend();
            matplot::save(path + "/plots/solutions" + ".png");
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