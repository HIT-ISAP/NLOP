#ifndef CONJUATEGRADIENTPARAMS_HPP
#define CONJUATEGRADIENTPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::ConjuateGradientParams
/// @brief Params of conjuate gradient optimization methods
class ConjuateGradientParams: public LineSearchParams
{
public:
    /// @brief Construct
    ConjuateGradientParams() { setDefaults(); }

    /// @brief Use default optimizer params
    void setDefaults()
    {
        max_iteration_times = 1000;
        iteration_times = 0;
        verbosity = SUMMARY;
        log_file = false;
        stepsize_method = WOLFEPOWELL;
        stepsize_params = new WolfePowellParams;
        min_gradient = 0.01;
    }

    /// @brief print params of optimizer
    void print() override
    {
        std::cout << "Conjuate Gradient Optimization" << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "stepsize method: " << StepsizeMethodTranslator(stepsize_method) << "\n";
        std::cout << "gradient thresthold: " << min_gradient << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Set and get mininum gradient to control iterations to stop
    void setMinGradient(const double value){ min_gradient = value; }
    double getMinGradient() const { return min_gradient; }

private:
    double min_gradient; // Gradient threshold to stop the iterations
};
}

#endif
