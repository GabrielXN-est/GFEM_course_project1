#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTIONS_H

#include "lin_alg.h"
#include <cmath>

class shape_functions: public Vector
{
    protected:
    double eval;
    int polinomial_order {};

    public:
    shape_functions(int po) : polinomial_order {po}, Vector{static_cast<std::size_t>(po+1)} {}

    void operator() (double e) {eval = e;}

    virtual double get_val (int index) =0;
    
    void mont_vector ()
    {
        for (int i {0}; i<=polinomial_order; i++)
        {
            vec[static_cast<size_t>(i)] = get_val(i);
        }
    }
};

class Mshape_functions_lag : public shape_functions
{
    public:
    double dist {2./static_cast<double>(polinomial_order)}; // distância entre os nós no elemento mestre

    Mshape_functions_lag(int po) : shape_functions {po}{}

    double get_val (int index);
};

class MDshape_functions_lag : public shape_functions
{
    public:
    double dist {2./static_cast<double>(polinomial_order)}; // distância entre os nós no elemento mestre

    MDshape_functions_lag(int po) : shape_functions {po}{}

    double get_val (int index);
};

class Mshape_functions_p_hier : public shape_functions
{
    public:

    Mshape_functions_p_hier(int po) : shape_functions {po}{}

    double get_val (int index);
};

class MDshape_functions_p_hier : public shape_functions
{
    public:

    MDshape_functions_p_hier(int po) : shape_functions {po}{}

    double get_val (int index);
};

double legendre_polynomy(double eval, int i);

double legendre_polynomy_integral(double eval, int i);

#endif