#ifndef RMSPROPPARAMS_HPP
#define RMSPROPPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::RMSPropParams
/// @brief Params of RMSProp optimization methods
class RMSPropParams: public LineSearchParams
{
public:
    /// @brief Constructor
    RMSPropParams() { setDefaults(); }

    /// @brief Set and get mininum gradient to control iterations to stop
    void setMinGradient(const double value) { min_gradient = value; }
    double getMinGradient() const { return min_gradient; }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        min_gradient = 0.01;
        alpha = 0.001;
        gamma = 0.9;
        epsilon = 1e-8;

        max_iteration_times = 10000;
        iteration_times = 0;
        verbosity = SUMMARY;
        log_file = false;
    }

    /// @brief Print params of optimizer
    void print() override
    {
        std::cout << "RMSProp Optimization" << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "gradient thresthold: " << min_gradient << "\n";
        std::cout << "learning rate: " << alpha << "\n";
        std::cout << "decay rate: " << gamma << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Setters and Getters of Params
    void setAlpha (const double value) { alpha = value; }
    void setGamma (const double value) { gamma = value; }
    void setEpsilon (const double value) { epsilon = value; }

    double getAlpha() const { return alpha; }
    double getGamma() const { return gamma; }
    double getEpsilon() const { return epsilon; }

private:
    double min_gradient; // gradient threshold to stop the iterations
    double alpha;        // learning rate
    double gamma;        // decay rate
    double epsilon;      // prevent division by 0
};
}

#endif
