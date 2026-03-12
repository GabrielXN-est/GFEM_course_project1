#include "mesh.h"
#include <iostream>
#include <iomanip>

void Mesh::debug_print_matrix_info()
{
    std::cout << "\n===== DEBUG: K_global_pos (antes eliminação) =====" << std::endl;
    std::cout << "Tamanho: " << K_global_pos.mat.size() << " x " << K_global_pos.mat[0].size() << std::endl;
    
    double min_val = 1e20, max_val = -1e20;
    int zeros = 0, nonzeros = 0;
    
    for (std::size_t i = 0; i < K_global_pos.mat.size(); i++)
    {
        for (std::size_t j = 0; j < K_global_pos.mat[i].size(); j++)
        {
            double val = K_global_pos[i][j];
            if (std::abs(val) < 1e-15) zeros++;
            else nonzeros++;
            
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
        }
    }
    
    std::cout << "Valores não-nulos: " << nonzeros << std::endl;
    std::cout << "Valores nulos: " << zeros << std::endl;
    std::cout << "Min: " << std::scientific << min_val << std::endl;
    std::cout << "Max: " << std::scientific << max_val << std::endl;
    
    std::cout << "\nDiagonal principal:" << std::endl;
    for (std::size_t i = 0; i < std::min(size_t(10), K_global_pos.mat.size()); i++)
    {
        std::cout << "K[" << i << "][" << i << "] = " << std::scientific << K_global_pos[i][i] << std::endl;
    }
    if (K_global_pos.mat.size() > 10)
    {
        for (std::size_t i = K_global_pos.mat.size() - 3; i < K_global_pos.mat.size(); i++)
        {
            std::cout << "K[" << i << "][" << i << "] = " << std::scientific << K_global_pos[i][i] << std::endl;
        }
    }
    
    std::cout << "\nElementos por linha (para detectar linhas nulas):" << std::endl;
    for (std::size_t i = 0; i < K_global_pos.mat.size(); i++)
    {
        int row_nonzeros = 0;
        for (std::size_t j = 0; j < K_global_pos.mat[i].size(); j++)
        {
            if (std::abs(K_global_pos[i][j]) > 1e-15) row_nonzeros++;
        }
        if (row_nonzeros == 0 || i < 3 || i >= K_global_pos.mat.size() - 3)
        {
            std::cout << "Linha " << i << ": " << row_nonzeros << " elementos não-nulos" << std::endl;
        }
    }
    
    try {
        double det = K_global_pos.determinant();
        std::cout << "\nDeterminante de K_global_pos: " << std::scientific << det << std::endl;
    } catch (const std::exception& e) {
        std::cout << "\nErro ao calcular determinante: " << e.what() << std::endl;
    }
}

