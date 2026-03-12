#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <type_traits>
#include <fstream>

#include "mesh.h"
#include "read_input.h"

int main()
{
    try {
    Mesh mesh {};
    read_input("/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_phier copy 2.txt", mesh);

    mesh.assemble_direct();
    mesh.solve();

    mesh.complete_U();
    std::cout << "Energy norm: " << mesh.energy_norm() << std::endl;
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