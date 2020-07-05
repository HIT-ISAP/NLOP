#ifndef LINESEARCHPARAMS_HPP
#define LINESEARCHPARAMS_HPP

#include <OptimizerParams/OptimizerParamsBase.hpp>
#include <StepsizeSearch/Accurate/AccurateSearchBase.hpp>


namespace NLOP {
/// @class NLOP::LineSearchParams
/// @brief Abstract class for params for all line search optimization
class LineSearchParams: public OptimizerParamsBase
{
public:
    // Stepsize search method
    enum StepsizeMethod{
        // Accuate stepsize search methods
        GOLDENSECTION,
        DICHOTOMOUS,
        FIBONACCI,
        // Inaccurate stepsize search methods
        ARMIJO,
        GOLDSTEIN,
        WOLFEPOWELL
    };
    /// @brief StepsizeMethod translator: from StepsizeMethod to std::string
    std::string StepsizeMethodTranslator(StepsizeMethod value)
    {
        std::string s;
        switch (value) {
        case GOLDENSECTION:
            s = "GOLDENSECTION";
            break;
        case DICHOTOMOUS:
            s = "DICHOTOMOUS";
            break;
        case FIBONACCI:
            s = "FIBONACCI";
            break;
        case ARMIJO:
            s = "ARMIJO";
            break;
        case GOLDSTEIN:
            s = "GOLDSTEIN";
            break;
        case WOLFEPOWELL:
            s = "WOLFEPOWELL";
            break;
        default:
            s = "UNDEFINED";
            break;
        }
        return s;
    }
    /// @brief StepsizeMethod translator: from std::string to StepsizeMethod
    StepsizeMethod StepsizeMethodTranslator(const std::string& value)
    {
        std::string s = value;
        if (s == "GOLDENSECTION")
            return GOLDENSECTION;
        if (s == "DICHOTOMOUS")
            return DICHOTOMOUS;
        if (s == "FIBONACCI")
            return FIBONACCI;
        if (s == "ARMIJO")
            return ARMIJO;
        if (s == "GOLDSTEIN")
            return GOLDSTEIN;
        if (s == "WOLFEPOWELL")
            return WOLFEPOWELL;

        return WOLFEPOWELL; // If invalid value is received, choose GOLDENSECTION mode
    }

    /// @brief Get and set stepsize method
    void setStepsizeMethod(std::string value) { stepsize_method = StepsizeMethodTranslator(value); }
    StepsizeMethod getStepsizeMethod() { return stepsize_method; }

protected:
    StepsizeMethod stepsize_method; // Stepsize search method (default GOLDENSECTION)

    virtual ~LineSearchParams() {}
};
}

#endif
