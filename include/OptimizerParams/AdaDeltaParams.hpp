#ifndef ADADELTAPARAMS_HPP
#define ADADELTAPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::AdaDeltaParams
/// @brief Params of AdaDelta optimization methods
class AdaDeltaParams: public LineSearchParams
{
public:
    /// @brief Constructor
    AdaDeltaParams() {}

    /// @brief Set mininum gradient to control iterations to stop
    void setMinGradient(const double value)
    {
        min_gradient = value;
    }

    double min_gradient = 0.01; // Gradient threshold to stop the iterations

    double alpha = 0.01; // learning rate
    double gamma = 0.9; // decay rate
    double epsilon = 1e-8; // Prevent division by 0
};
}


#endif
