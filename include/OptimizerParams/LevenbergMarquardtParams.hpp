#ifndef LM_PARAMS_HPP
#define LM_PARAMS_HPP

#include <OptimizerParams/TrustRegionParams.hpp>

namespace NLOP {
/// @class NLOP::LevenbergMarquardtParams
/// @brief Params of Levenberg-Marquardt optimization methods
class LevenbergMarquardtParams: public LineSearchParams
{
public:
    /// @brief Constructor
    LevenbergMarquardtParams() {}

    double min_delta_x = 0.0001; // Gradient threshold to stop the iterations
    double init_epsilon = 100;
};
}

#endif
