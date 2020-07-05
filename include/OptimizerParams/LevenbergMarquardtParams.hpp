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
    LevenbergMarquardtParams() { setDefaults(); }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        max_iteration_times = 1000;
        iteration_times = 0;
        verbosity = SUMMARY;
        min_delta_x = 0.0001;
        init_epsilon = 4;
    }

    /// @brief print params of optimizer
    void print(const std::string &str) override
    {
        std::cout << str << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "stepsize thresthold: " << min_delta_x << "\n";
        std::cout << "initial damping factor: " << init_epsilon << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Setters and Getters of params
    void setMinDeltaX(const double value) { min_delta_x = value; }
    void setInitEpsilon(const double value) { init_epsilon = value; }

    double getMinDeltaX() const { return min_delta_x; }
    double getInitEpsilon() const { return init_epsilon; }

private:
    double min_delta_x = 0.0001; // Gradient threshold to stop the iterations
    double init_epsilon = 4;     // Initial damping factor
};
}

#endif
