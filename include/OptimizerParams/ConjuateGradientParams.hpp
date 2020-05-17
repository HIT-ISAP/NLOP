#ifndef CONJUATEGRADIENTPARAMS_HPP
#define CONJUATEGRADIENTPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::ConjuateGradientParams
/// @brief Params of conjuate gradient optimization methods
class ConjuateGradientParams: public LineSearchParams
{
public:
    /// @brief Construct
    ConjuateGradientParams() {}

    /// @brief Use default optimizer params
    void setDefaults()
    {
        max_iteration_times = 10000;
        iteration_times = 0;
        verbosity = false;
        stepsize_search_method = GOLDENSECTION;
        min_gradient = 0.01;
    }

    /// @brief Set mininum gradient to control iterations to stop
    void setMinGradient(const double value)
    {
        min_gradient = value;
    }

    double min_gradient = 0.01; // Gradient threshold to stop the iterations
};
}

#endif
