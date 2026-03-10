#include <cmath>
#include <stdexcept>

#include "integration_points.h"

integration_points Gauss_quad_points (int order)
{
    int points {static_cast<int>(ceil((static_cast<double>(order)+1)/2.))};
    switch (points)
    {
        case 1:
            return integration_points{{0.0}, {2.0}};
        case 2:
            return integration_points{{1./sqrt(3.), -1./sqrt(3.)}, {1.0, 1.0}};
        case 3:
            return integration_points{{0., -sqrt(3./5.), sqrt(3./5.)}, {8./9., 5./9., 5./9.}};
        case 4:
            return integration_points{{-sqrt(3./7.-2./7.*sqrt(6./5.)), sqrt(3./7.-2./7.*sqrt(6./5.)), 
                                       -sqrt(3./7.+2./7.*sqrt(6./5.)), sqrt(3./7.+2./7.*sqrt(6./5.))},
                                        {(18.+sqrt(30.))/36., (18.+sqrt(30.))/36., (18.-sqrt(30.))/36., (18.-sqrt(30.))/36.}};
        case 5:
            return integration_points{{0., -1./3.*sqrt(5.-2.*sqrt(10./7.)), 1./3.*sqrt(5.-2.*sqrt(10./7.)), 
                                       -1./3.*sqrt(5.+2.*sqrt(10./7.)), 1./3.*sqrt(5.+2.*sqrt(10./7.))},
                                        {128./225., (322.+13.*sqrt(70.))/900., (322.+13.*sqrt(70.))/900., (322.-13.*sqrt(70.))/900., (322.-13.*sqrt(70.))/900.}};
        case 6:
            return integration_points{{-0.9324695142, 0.9324695142, 0.6612093865, -0.6612093865, 0.2386191861, -0.2386191861},
                                        {0.1713244924, 0.1713244924, 0.3607615730, 0.3607615730, 0.4679139346, 0.4679139346}};
        case 7:
            return integration_points{{0., -0.9491079123, 0.9491079123, -0.7415311856, 0.7415311856, -0.4058451514, 0.4058451514},
                                        {0.4179591837, 0.1294849662, 0.1294849662, 0.2797053915, 0.2797053915, 0.3818300505, 0.3818300505}};
        case 8:
            return integration_points{{0.1834346425, -0.1834346425,0.5255324099, -0.5255324099,
                                        0.7966664774, -0.7966664774, 0.9602898565, -0.9602898565},
                                        {0.3626837834,0.3626837834,0.3137066459, 0.3137066459,
                                        0.2223810345, 0.2223810345, 0.1012285363, 0.1012285363}};
        default:
            throw std::invalid_argument("Unsupported order for Gauss quadrature");
    }
}