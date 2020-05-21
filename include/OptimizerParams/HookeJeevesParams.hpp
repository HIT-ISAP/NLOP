#ifndef HOOKEJEEVESPARAMS_HPP
#define HOOKEJEEVESPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::HookeJeevesParams
/// @brief Params of Hooke Jeeves optimization methods
class HookeJeevesParams: public LineSearchParams
{
public:
    /// @brief Constructor
    HookeJeevesParams() {}

    double alpha = 1;
    double init_delta = 0.2;
    double epsilon = 0.001;
};
}

#endif
