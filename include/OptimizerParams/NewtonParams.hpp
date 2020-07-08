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
    NewtonParams() { setDefaults(); }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        max_iteration_times = 1000;
        iteration_times = 0;
        verbosity = SUMMARY;
        log_file = false;
        min_delta_x = 0.01;
    }

    /// @brief print params of optimizer
    void print() override
    {
        std::cout << "Newton's Method Optimization" << std::endl;
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "stepsize thresthold: " << min_delta_x << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Setters and Getters of params
    void setMinDeltaX(const double value) { min_delta_x = value; }
    double getMinDeltaX() const { return min_delta_x; }

private:
    double min_delta_x; // stepsize threshold to stop the iterations
};
}

#endif
