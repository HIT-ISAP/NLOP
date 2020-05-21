#ifndef ADAGRADPARAMS_HPP
#define ADAGRADPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::AdagradParams
/// @brief Params of Adagrad optimization methods
class AdagradParams: public LineSearchParams
{
public:
    /// @brief Constructor
    AdagradParams() {}

    /// @brief Set mininum gradient to control iterations to stop
    void setMinGradient(const double value)
    {
        min_gradient = value;
    }

    double min_gradient = 0.01; // Gradient threshold to stop the iterations

    double alpha = 0.1; // learning rate
    double epsilon = 1e-8; // Prevent division by 0
};
}

#endif
