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
    read_input("/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input.txt", mesh);

    mesh.assemble_penalty();
    mesh.solve();

    std::cout << "Energy norm: " << mesh.energy_norm() << std::endl;
    return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}