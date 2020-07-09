#ifndef INACCURATE_SEARCH_PARAMS_HPP
#define INACCURATE_SEARCH_PARAMS_HPP

#include <StepsizeSearchParams/StepsizeSearchParamsBase.hpp>

namespace NLOP {
/// @class NLOP::InaccurateSearchParams
/// @brief Abstract base class for params of all accurate stepsizes search method
class InaccurateSearchParams: public StepsizeSearchParamsBase
{
public:
    /// @brief Constructor
    InaccurateSearchParams() { setDefaults(); }

    /// @brief Use default stepsize searcher params
    void setDefaults() override
    {
        alpha = 1.5;
        beta = 0.5;
        init_lambda_factor = 0.1;
    }

    /// @brief Getters and Setters for params
    void setIncreaseFactor(const double value) override { alpha = value; }
    void setDecreaseFactor(const double value) override { beta = value; }
    void setInitLambdaFactor(const double value) override { init_lambda_factor = value; }

    double getIncreaseFactor() const override { return alpha; }
    double getDecreaseFactor() const override { return beta; }
    double getInitLambdaFactor() const override { return init_lambda_factor; }

protected:
    double alpha;               // increase factor
    double beta;                // decrease factor
    double init_lambda_factor;  // initial stepsize factor

};
}

#endif
