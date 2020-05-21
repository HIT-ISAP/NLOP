#ifndef ADAMPARAMS_HPP
#define ADAMPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::AdamParams
/// @brief Params of Adam optimization methods
class AdamParams: public LineSearchParams
{
public:
    /// @brief Constructor
    AdamParams() {}

    /// @brief Set mininum gradient to control iterations to stop
    void setMinGradient(const double value)
    {
        min_gradient = value;
    }

    double min_gradient = 0.01; // Gradient threshold to stop the iterations

    double alpha = 0.01; // learning rate
    double gamma_v = 0.9; // decay rate v
    double gamma_s = 0.999; // decay rate s
    double epsilon = 1e-8; // Prevent division by 0
};
}

#endif
