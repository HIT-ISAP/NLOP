#ifndef ADAMPARAMS_HPP
#define ADAMPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::AdamParams
/// @brief Params of Adam optimization methods
class AdamParams: public LineSearchParams
{
public:
    /// @brief Constructor
    AdamParams() { setDefaults(); }

    /// @brief Set and get mininum gradient to control iterations to stop
    void setMinGradient(const double value) { min_gradient = value; }
    double getMinGradient() const { return min_gradient; }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        min_gradient = 0.01;
        alpha = 0.01;
        gamma_v = 0.9;
        gamma_s = 0.999;
        epsilon = 1e-8;

        max_iteration_times = 10000;
        iteration_times = 0;
        verbosity = SUMMARY;
        log_file = false;
    }

    /// @brief print params of optimizer
    void print() override
    {
        std::cout << "Adam Optimization" << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "gradient thresthold: " << min_gradient << "\n";
        std::cout << "learning rate: " << alpha << "\n";
        std::cout << "decay rate v: " << gamma_v << "\n";
        std::cout << "decay rate s: " << gamma_s << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Setters and Getters of params
    void setAlpha(const double value) { alpha = value; }
    void setGammaV(const double value) { gamma_v = value; }
    void setGammaS(const double value) { gamma_s = value; }
    void setEpsilon(const double value) { epsilon = value; }

    double getAlpha() const { return alpha; }
    double getGammaV() const { return gamma_v; }
    double getGammaS() const { return gamma_s; }
    double getEpsilon() const { return epsilon; }

private:
    double min_gradient = 0.01; // Gradient threshold to stop the iterations
    double alpha = 0.01;        // Learning rate
    double gamma_v = 0.9;       // Decay rate v
    double gamma_s = 0.999;     // Decay rate s
    double epsilon = 1e-8;      // Prevent division by 0
};
}

#endif
