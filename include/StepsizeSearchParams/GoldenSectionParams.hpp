#ifndef GOLDENSECTION_PARAMS_HPP
#define GOLDENSECTION_PARAMS_HPP

#include <StepsizeSearchParams/AccurateSearchParams.hpp>

namespace NLOP {
/// @class NLOP::GoldenSectionParams
/// @brief params of golden section method
class GoldenSectionParams: public AccurateSearchParams
{
public:
    /// @brief Constructor

    /// @brief Use default stepsize searcher params
    void setDefaults() override
    {
        setLowerBound(0);
        setUpperBound(1);
        setMaxIterations(20);
        setEpsilon(0.001);
    }
};
}

#endif
