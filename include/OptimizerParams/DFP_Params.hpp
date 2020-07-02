#ifndef DFP_PARAMS_HPP
#define DFP_PARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::DFP_Params
/// @brief Params of DFP Quasi-Newton method
class DFP_Params: public LineSearchParams
{
public:
    /// @brief Constructor
    DFP_Params() {}

    double min_gradient = 0.01; // gradient threshold to stop iterations
};
}

#endif
