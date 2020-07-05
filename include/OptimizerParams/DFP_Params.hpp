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
    DFP_Params() { setDefaults(); }

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
