#ifndef BFGS_PARAMS_HPP
#define BFGS_PARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::BFGS_Params
/// @brief Params of BFGS Quasi-Newton method
class BFGS_Params: public LineSearchParams
{
public:
    /// @brief Constructor
    BFGS_Params() {}

    double min_gradient = 0.01; // gradient threshold to stop iterations
};
}

#endif
