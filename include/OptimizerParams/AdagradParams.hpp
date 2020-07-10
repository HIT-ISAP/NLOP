#ifndef ADAGRADPARAMS_HPP
#define ADAGRADPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::AdagradParams
/// @brief Params of Adagrad optimization methods
class AdagradParams: public LineSearchParams
{
public:
    /// @brief Constructor
    AdagradParams() { setDefaults(); }

    /// @brief Set and get mininum gradient to control iterations to stop
    void setMinGradient(const double value) { min_gradient = value; }
    double getMinGradient() const { return min_gradient; }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        min_gradient = 0.01;
        alpha = 0.1;
        epsilon = 1e-8;

        max_iteration_times = 10000;
        iteration_times = 0;
        verbosity = SUMMARY;
        log_file = false;
    }

    /// @brief Print params of optimizer
    void print() override
    {
        std::cout << "Adagrad Optimization" << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "gradient thresthold: " << min_gradient << "\n";
        std::cout << "learning rate: " << alpha << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Setters and Getters of params
    void setAlpha(const double value) { alpha = value; }
    void setEpsilon(const double value) { epsilon = value; }

    double getAlpha() const { return alpha; }
    double getEpsilon() const { return epsilon; }

private:
    double min_gradient; // gradient threshold to stop the iterations
    double alpha;        // learning rate
    double epsilon;      // prevent division by 0
};
}

#endif
