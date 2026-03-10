#ifndef INTEGRATION_POINTS_H
#define INTEGRATION_POINTS_H

#include <vector>

struct integration_points
{
    std::vector<double> points {};
    std::vector<double> weights {};
};

integration_points Gauss_quad_points (int order);

#endif