#ifndef NESTEROVMOMENTUMPARAMS_HPP
#define NESTEROVMOMENTUMPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::NesterovMomentumParams
/// @brief Params of Momentum optimization methods
class NesterovMomentumParams: public LineSearchParams
{
public:
    /// @brief Constructor
    NesterovMomentumParams() { setDefaults(); }

    /// @brief Set and get mininum gradient to control iterations to stop
    void setMinGradient(const double value) { min_gradient = value; }
    double getMinGradient() const { return min_gradient; }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        min_gradient = 0.01;
        alpha = 0.001;
        beta = 0.9;

        max_iteration_times = 10000;
        iteration_times = 0;
        verbosity = SUMMARY;
        log_file = false;
    }

    /// @brief print params of optimizer
    void print() override
    {
        std::cout << "Nesterov Momentum Optimization" << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "gradient thresthold: " << min_gradient << "\n";
        std::cout << "Learning rate: " << alpha << "\n";
        std::cout << "momentum hyper-parameter: " << beta << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Setters and Getters of params
    void setAlpha(const double value) { alpha = value; }
    void setBeta(const double value) { beta = value; }

    double getAlpha() const { return alpha; }
    double getBeta() const { return beta; }

private:
    double min_gradient;  // Gradient threshold to stop the iterations
    double alpha;         // Learning rate
    double beta;          // Momentum hyper-parameter

};
}

#endif
