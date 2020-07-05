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
    BFGS_Params() { setDefaults(); }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        max_iteration_times = 1000;
        iteration_times = 0;
        verbosity = SUMMARY;
        stepsize_method = WOLFEPOWELL;
        min_gradient = 0.01;
    }

    /// @brief print params of optimizer
    void print(const std::string &str) override
    {
        std::cout << str << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "stepsize method: " << StepsizeMethodTranslator(stepsize_method) << "\n";
        std::cout << "gradient thresthold: " << min_gradient << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Set and get mininum gradient to control iterations to stop
    void setMinGradient(const double value) { min_gradient = value; }
    double getMinGradient() const { return min_gradient; }

private:
    double min_gradient = 0.01; // gradient threshold to stop iterations
};
}

#endif
