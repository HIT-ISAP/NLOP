#ifndef ADADELTAPARAMS_HPP
#define ADADELTAPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::AdaDeltaParams
/// @brief Params of AdaDelta optimization methods
class AdaDeltaParams: public LineSearchParams
{
public:
    /// @brief Constructor
    AdaDeltaParams() { setDefaults(); }

    /// @brief Set and get mininum gradient to control iterations to stop
    void setMinGradient(const double value) { min_gradient = value; }
    double getMinGradient() const { return min_gradient; }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        min_gradient = 0.01;
        alpha = 0.01;
        gamma = 0.9;
        epsilon = 1e-8;
        sgd_times = 1000;
        stepsize_method = WOLFEPOWELL;
        stepsize_params = new WolfePowellParams;

        max_iteration_times = 10000;
        iteration_times = 0;
        verbosity = SUMMARY;
        log_file = false;
    }

    /// @brief print params of optimizer
    void print() override
    {
        std::cout << "AdaDelta Optimization" << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "gradient thresthold: " << min_gradient << "\n";
        std::cout << "SGD iteration times: " << sgd_times << "\n";
        std::cout << "learning rate: " << alpha << "\n";
        std::cout << "decay rate: " << gamma << "\n";
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Setters and Getters of Params
    void setAlpha (const double value) { alpha = value; }
    void setGamma (const double value) { gamma = value; }
    void setEpsilon (const double value) { epsilon = value; }
    void setSGDTimes(const size_t value) { sgd_times = value; }

    double getAlpha() const { return alpha; }
    double getGamma() const { return gamma; }
    double getEpsilon() const { return epsilon; }
    size_t getSGDTimes() const { return sgd_times; }

private:
    size_t sgd_times;    // Perform several times of steepest gradient descent to launch the AdaDelta method (default = 1000)
    double min_gradient; // Gradient threshold to stop the iterations (default = 0.01)
    double alpha;        // Learning rate (default = 0.01)
    double gamma;        // Decay rate (default = 0.9)
    double epsilon;      // Prevent division by 0 (default = 1e-8)
};
}


#endif
