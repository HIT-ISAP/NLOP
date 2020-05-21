#ifndef MOMENTUMPARAMS_HPP
#define MOMENTUMPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::MomentumParams
/// @brief Params of Momentum optimization methods
class MomentumParams: public LineSearchParams
{
public:
    /// @brief Constructor
    MomentumParams() {}

    /// @brief Set mininum gradient to control iterations to stop
    void setMinGradient(const double value)
    {
        min_gradient = value;
    }

    double min_gradient = 0.01; // Gradient threshold to stop the iterations

    double alpha = 0.003; // learning rate
    double beta = 0.9; // momentum hyper-parameters
};
}

#endif
