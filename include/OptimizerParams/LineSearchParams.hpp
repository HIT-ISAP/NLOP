#ifndef LINESEARCHPARAMS_HPP
#define LINESEARCHPARAMS_HPP

#include <OptimizerParams/OptimizerParamsBase.hpp>
#include <StepsizeSearch/Accurate/AccurateSearchBase.hpp>


namespace NLOP {

/// @class NLOP::LineSearchParams
/// @brief Abstract class for params for all line search optimization
/// @param VariableType The type of variable x

class LineSearchParams: public OptimizerParamsBase
{
public:
    /// all one dimensional search options
    enum OneDimensionalSearchMethod {
        BISECTION,
        GOLDENSECTION,
        FIBONACCI,
        DICHOTOMOUS,
        NEWTON
    };

    /// @brief Set one dimensional search method
    void setOneDimensionalSearchMethod(const std::string& searcher)
    {
        if (searcher == "BISECTION")
            stepsize_search_method = BISECTION;
        else if (searcher == "GOLDENSECTION")
            stepsize_search_method = GOLDENSECTION;
        else if (searcher == "FIBONACCI")
            stepsize_search_method = FIBONACCI;
        else if (searcher == "DICHOTOMOUS")
            stepsize_search_method = DICHOTOMOUS;
        else if (searcher == "NEWTON")
            stepsize_search_method = NEWTON;
        else
            std::cerr << "Invalid Option" << std::endl;
    }
protected:
    virtual ~LineSearchParams() {}

    OneDimensionalSearchMethod stepsize_search_method = GOLDENSECTION;


};
}

#endif
