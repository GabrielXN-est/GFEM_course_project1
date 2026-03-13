#include "shape_functions.h"


double Mshape_functions_lag::get_val (int index)
{
    double Ni {1};
    for (int r {0}; r<=polinomial_order; r++)
    {
        if (r != index)
            {Ni *= (eval-dist*r+1)/(dist*index-dist*r);}
    }
    return Ni;
}

double MDshape_functions_lag::get_val (int index)
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

double Mshape_functions_p_hier::get_val (int index) // index 0-based no vetor
{
    if (index == 0)
        return (1.-eval)/2.;
    else if (index == polinomial_order)
        return (1.+eval)/2.;
    else
        return std::sqrt((2.* static_cast<double>(index+2) -3)/2.)*(legendre_polynomy_integral(eval, index)-legendre_polynomy_integral(-1, index));
}

double MDshape_functions_p_hier::get_val (int index)
{
    if (index == 0)
        return -1./2.;
    else if (index == polinomial_order)
        return 1./2.;
    else
        return std::sqrt((2.* static_cast<double>(index+2) -3)/2.)*(legendre_polynomy(eval, index));
}

double legendre_polynomy(double eval, int i)
{
    switch (i)
    {
    case 0:
        return 1.;
    case 1:
        return eval;
    case 2:
        return (3.*eval*eval-1.)/2.;
    case 3:
        return (-3.+5.*std::pow(eval, 2))*eval/2.;
    case 4:
        return (3.-30.*std::pow(eval, 2)+35.*std::pow(eval, 4))/8.;
    case 5:
        return (15.-70.*std::pow(eval, 2)+63.*std::pow(eval, 4))*eval/8.;
    case 6:
        return (-5.+105.*std::pow(eval, 2)-315.*std::pow(eval, 4)+231.*std::pow(eval, 6))/16.;
    case 7:
        return (-35.+315.*std::pow(eval, 2)-693.*std::pow(eval, 4)+429.*std::pow(eval, 6))*eval/16.;
    case 8:
        return (35.-1260.*std::pow(eval, 2)+6930.*std::pow(eval, 4)-12012.*std::pow(eval, 6)+6435.*std::pow(eval, 8))/128.;
    default:
        throw std::invalid_argument("Legendre polynomial order not implemented");
    }
}

double legendre_polynomy_integral(double eval, int i)
{
    switch (i)
    {
    case 0:
        return eval;
    case 1:
        return (eval*eval)/2.;
    case 2:
        return (-eval+std::pow(eval, 3))/2.;
    case 3:
        return (std::pow(eval, 2)*(-6.+5.*std::pow(eval, 2)))/8.;
    case 4:
        return eval*(3.-10.*std::pow(eval, 2)+7.*std::pow(eval, 4))/8.;
    case 5:
        return std::pow(eval, 2)*(5.-35.*std::pow(eval, 2)+21.*std::pow(eval, 4))/16.;
    case 6:
        return eval*(-5.+35. *std::pow(eval, 2)-63.*std::pow(eval, 4)+33.*std::pow(eval, 6))/16.;
    case 7:
        return std::pow(eval, 2)*(-140.+630.*std::pow(eval, 2)-924.*std::pow(eval, 4)+429.*std::pow(eval, 6))/128.;
    case 8:
        return eval*(35.-420.*std::pow(eval, 2)+1386.*std::pow(eval, 4)-1716.*std::pow(eval, 6)+715.*std::pow(eval, 8))/128.;
    default:
        throw std::invalid_argument("Legendre polynomial order not implemented");
    }
}