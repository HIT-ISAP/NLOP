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
    ConjuateGradientParams() {}

    // Approximate method to avoid computing Hessian
    enum ApproximateMethod {
        FR, // Fletcher-Reeves
        PR // Polak-Ribiere
    };

    /// @brief Use default optimizer params
    void setDefaults()
    {
        max_iteration_times = 10000;
        iteration_times = 0;
        verbosity = false;
        stepsize_search_method = GOLDENSECTION;
        min_gradient = 0.01;
        approximation = FR;
    }

    void setApproximateMethod(std::string method)
    {
        approximation = enumTranslator(method);
    }

    ApproximateMethod enumTranslator(std::string method)
    {
        ApproximateMethod result;
        if (method == "FR")
            result = FR;
        else if (method == "PR")
            result = PR;
        else
        {
            std::cerr << "Values ​​other than FR and PR are invalid!" << std::endl;
            result = FR;
        }
        return result;
    }

    /// @brief Set mininum gradient to control iterations to stop
    void setMinGradient(const double value)
    {
        min_gradient = value;
    }

    double min_gradient = 0.01; // Gradient threshold to stop the iterations
    ApproximateMethod approximation = FR; // Approximate method to avoid computing Hessian
                                          // Default == FR

};
}

#endif
