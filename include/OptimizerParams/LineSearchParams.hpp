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
        // Accurate stepsize search methods
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

    /// @brief Getters and setters for common params
    void setStepsizeMethod(std::string value) { stepsize_method = StepsizeMethodTranslator(value); }
    void setStepsizeUpperBound(const double value) { stepsize_params->setUpperBound(value); }
    void setStepsizeLowerBound(const double value) { stepsize_params->setLowerBound(value); }
    void setStepsizeMaxIterations(const size_t value) { stepsize_params->setMaxIterations(value); }

    StepsizeMethod getStepsizeMethod() { return stepsize_method; }
    double getStepsizeUpperBound() const { return stepsize_params->getUpperBound(); }
    double getStepsizeLowerBound() const { return stepsize_params->getLowerBound(); }
    size_t getStepsizeMaxIterations() const { return stepsize_params->getMaxIterations(); }

    /// @brief Getters and Setters for params specialized for accurate searching method
    void setStepsizeAccuracy(const double value)
    {
        if (stepsize_method == GOLDENSECTION || stepsize_method == FIBONACCI || stepsize_method == DICHOTOMOUS)
            stepsize_params->setStepsizeAccuracy(value);
        else
            std::cerr << "The selected stepsize searching method is an inaccurate type, attribute Stepsize Accuracy is not existing" << std::endl;
    }
    double getStepsizeAccuracy() const
    {
        if (stepsize_method == GOLDENSECTION || stepsize_method == FIBONACCI || stepsize_method == DICHOTOMOUS)
            return stepsize_params->getStepsizeAccuracy();
        else
            std::cerr << "The selected stepsize searching method is an inaccurate type, attribute Stepsize Accuracy is not existing" << std::endl;
    }

    /// @brief Getters and Setters for params specialized for inaccurate searching method
    void setStepsizeIncreaseFactor(const double value)
    {
        if (stepsize_method == GOLDSTEIN || stepsize_method == ARMIJO || stepsize_method == WOLFEPOWELL)
            stepsize_params->setIncreaseFactor(value);
        else
            std::cerr << "The selected stepsize searching method is an accurate type, attribute Increase Factor is not existing" << std::endl;
    }
    void setStepsizeDecreaseFactor(const double value)
    {
        if (stepsize_method == GOLDSTEIN || stepsize_method == ARMIJO || stepsize_method == WOLFEPOWELL)
            stepsize_params->setDecreaseFactor(value);
        else
            std::cerr << "The selected stepsize searching method is an accurate type, attribute Increase Factor is not existing" << std::endl;
    }
    double getStepsizeIncreaseFactor() const
    {
        if (stepsize_method == GOLDSTEIN || stepsize_method == ARMIJO || stepsize_method == WOLFEPOWELL)
            return stepsize_params->getIncreaseFactor();
        else
            std::cerr << "The selected stepsize searching method is an accurate type, attribute Dncrease Factor is not existing" << std::endl;
    }
    double getStepsizeDecreaseFactor() const
    {
        if (stepsize_method == GOLDSTEIN || stepsize_method == ARMIJO || stepsize_method == WOLFEPOWELL)
            return stepsize_params->getDecreaseFactor();
        else
            std::cerr << "The selected stepsize searching method is an accurate type, attribute Dncrease Factor is not existing" << std::endl;
    }

    StepsizeSearchParamsBase* stepsize_params;  // params for stepsize searching

protected:
    StepsizeMethod stepsize_method;             // stepsize searching method (default WOLFEPOWELL)


    virtual ~LineSearchParams() {}
};
}

#endif
