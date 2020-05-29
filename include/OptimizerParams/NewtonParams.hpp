#ifndef NEWTONPARAMS_HPP
#define NEWTONPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::NewtonParams
/// @brief Params of Newton optimization methods
class NewtonParams: public LineSearchParams
{
public:
    /// @brief Constructor
    NewtonParams() {}

    double min_delta_x = 0.01; // stepsize threshold to stop the iterations
};
}

#endif
