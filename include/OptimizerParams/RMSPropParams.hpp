#ifndef RMSPROPPARAMS_HPP
#define RMSPROPPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::RMSPropParams
/// @brief Params of RMSProp optimization methods
class RMSPropParams: public LineSearchParams
{
public:
    /// @brief Constructor
    RMSPropParams() {}

    /// @brief Set mininum gradient to control iterations to stop
    void setMinGradient(const double value)
    {
        min_gradient = value;
    }

    double min_gradient = 0.01; // Gradient threshold to stop the iterations

    double alpha = 0.001; // learning rate
    double gamma = 0.9; // decay rate
    double epsilon = 1e-8; // Prevent division by 0
};
}

#endif
