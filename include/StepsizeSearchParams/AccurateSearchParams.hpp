#ifndef ACCURATE_SEARCH_PARAMS_HPP
#define ACCURATE_SEARCH_PARAMS_HPP

#include <StepsizeSearchParams/StepsizeSearchParamsBase.hpp>

namespace NLOP {
/// @class NLOP::AccurateSearchParams
/// @brief Abstract base class for params of all accurate stepsizes search method
class AccurateSearchParams: public StepsizeSearchParamsBase
{
public:
    /// @brief Constructor
    AccurateSearchParams() { setDefaults(); }

    /// @brief Use default stepsize searcher params
    void setDefaults() override
    {
        setLowerBound(0);
        setUpperBound(1);
        setMaxIterations(20);
        setStepsizeAccuracy(0.001);
    }

    /// @brief Getters and Setters for params
    void setStepsizeAccuracy(const double value) override { epsilon = value; }
    double getStepsizeAccuracy() const override { return epsilon; }

protected:
    double epsilon = 0.001;     // stepsize accuracy
};
}

#endif
