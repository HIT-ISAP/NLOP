#ifndef HOOKEJEEVESPARAMS_HPP
#define HOOKEJEEVESPARAMS_HPP

#include <OptimizerParams/LineSearchParams.hpp>

namespace NLOP {
/// @class NLOP::HookeJeevesParams
/// @brief Params of Hooke Jeeves optimization methods
class HookeJeevesParams: public LineSearchParams
{
public:
    /// @brief Constructor
    HookeJeevesParams() { setDefaults(); }

    /// @brief Use default optimizer params
    void setDefaults() override
    {
        alpha = 1;
        init_delta = 0.2;
        epsilon = 0.001;

        max_iteration_times = 10000;
        iteration_times = 0;
        verbosity = SUMMARY;
        log_file = false;
    }

    /// @brief print params of optimizer
    void print() override
    {
        std::cout << "Hooke&Jeeves Optimization" << "\n";
        std::cout << "*********************************************" << "\n";
        std::cout << "maximum iterations: " << max_iteration_times << "\n";
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << "\n";
        std::cout << "initial stepsize: " << init_delta << "\n";
        std::cout << "acceleration factor: " << alpha << "\n";
        std::cout << "stepsize thresthold: " << epsilon << "\n";;
        std::cout << "*********************************************" << std::endl;
    }

    /// @brief Setters and Getters of params
    void setAlpha(const double value) { alpha = value; }
    void setInitStepsize(const double value) { init_delta = value; }
    void setEpsilon(const double value) { epsilon = value; }

    double getAlpha() const { return alpha; }
    double getInitStepsize() const { return init_delta; }
    double getEpsilon() const { return epsilon; }

private:
    double alpha = 1;         // Acceleration factor
    double init_delta = 0.2;  // Initial stepsize
    double epsilon = 0.001;   // Stepsize threshold to stop the iteration
};
}

#endif
