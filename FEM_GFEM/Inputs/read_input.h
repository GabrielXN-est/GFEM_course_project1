#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <vector>
#include <string>

#include "body_func.h"

class Mesh;

class Properties
{
    public:
    int id;
    double A {};
    std::vector<double> E {};
    std::vector<double> Exlim {};
    double C {};
    body_functions* bf_func;

    void clear()
    {
        id = 0;
        A = C = 0;
        E.clear();
        delete bf_func;
        bf_func = nullptr;
    }

    ~Properties() {delete bf_func;}
};

void read_input (const std::string& filename, Mesh& mesh);

#endif