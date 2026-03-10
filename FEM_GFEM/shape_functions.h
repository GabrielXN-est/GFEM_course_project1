#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTIONS_H

#include "lin_alg.h"

class shape_functions: public Vector
{
    protected:
    double eval;
    int polinomial_order {};

    public:
    double dist {2./static_cast<double>(polinomial_order)};
    shape_functions(int po) : polinomial_order {po}, Vector{static_cast<std::size_t>(po+1)} {}

    void operator() (double e)
    {
        eval = e;
    }

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

    Mshape_functions_lag(int po) : shape_functions {po}{}

    double get_val (int index)
    {
        double Ni {1};
        for (int r {0}; r<=polinomial_order; r++)
        {
            if (r != index)
                {Ni *= (eval-dist*r+1)/(dist*index-dist*r);}
        }
        return Ni;
    }
};

class MDshape_functions_lag : public shape_functions
{
    public:
    
    MDshape_functions_lag(int po) : shape_functions {po}{}

    double get_val (int index)
    {
        double Ni {0};
        double temp {};
        for (int r {0}; r<=polinomial_order; r++)
        {
            if (r != index)
                {
                    temp = 1/(dist*index-dist*r);
                    for (int j {0}; j<=polinomial_order; j++)
                    {
                        if (j != index && j != r)
                            {temp *= (eval-j*dist+1)/(dist*index-dist*j);}
                    }
                    Ni += temp;
                }
        }
        return Ni;
    }
};

#endif