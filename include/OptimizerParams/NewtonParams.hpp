#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::NewtonParams
/// @brief Params of Newton optimization methods
class NewtonParams
{
public:
    /// @brief Constructor
    NewtonParams() {}

    /// @brief Set mininum gradient to control iterations to stop
    void setMinGradient(const double value)
    {
        min_gradient = value;
    }

    double min_gradient = 0.01; // Gradient threshold to stop the iterations
};
}
