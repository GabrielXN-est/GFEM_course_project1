#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <type_traits>
#include <fstream>

#include "create_input.h"
#include "mesh.h"
#include "read_input.h"

int main()
{
    try {
    std::string filename {"/home/labmec/Downloads/GFEM Course/Projects/Projeto 1/input_generated.txt"};

    generate_input(filename, 2, 2, "lBar", 1, std::vector<double> {1}, std::vector<double> {}, 1, 0,
    std::vector<double> {{0.,0.}}, std::vector<int> {0, 1}, std::vector<double> {}, std::vector<int> {},
    10, 0.5, 0.2);

    Mesh mesh {};
    read_input(filename, mesh);

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