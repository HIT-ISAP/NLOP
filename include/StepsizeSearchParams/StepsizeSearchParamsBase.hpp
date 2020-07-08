#ifndef STEPSIZE_SEARCH_PARAMS_HPP
#define STEPSIZE_SEARCH_PARAMS_HPP

#include <string>
#include <iostream>

namespace NLOP {
/// @class NLOP::StepsizeSearchParamsBase
/// @brief Abstract base class for params of all stepsizes search methods
class StepsizeSearchParamsBase
{
public:
    /// @Constructor
    StepsizeSearchParamsBase() { setDefaults(); }

    /// @brief Use default stepsize searcher params
    virtual void setDefaults()
    {
        iteration_times = 0;
        max_iteration_times = 20;
        upper_bound = 1;
        lower_bound = 0;
    }

    /// @brief iteration_times += 1
    void nextIteration() { iteration_times++; }

    /// @brief Getters and Setters for common params
    size_t getMaxIterations() const { return max_iteration_times; }
    size_t getIterationTimes() const { return iteration_times; }
    double getUpperBound() const { return upper_bound; }
    double getLowerBound() const { return lower_bound; }

    void setMaxIterations(const size_t value) { max_iteration_times = value; }
    void setIterationTimes(const size_t value) { iteration_times = value; }
    void setUpperBound(const double value) { upper_bound = value; }
    void setLowerBound(const double value) { lower_bound = value; }

    /// @brief Getters and Setters for params specialized for accurate searching method
    virtual void setStepsizeAccuracy(const double value) {}
    virtual double getStepsizeAccuracy() {}

    /// @brief Getters and Setters for params specialized for inaccurate searching method
    virtual void setIncreaseFactor(const double value) {}
    virtual void setDecreaseFactor(const double value) {}
    virtual void setInitLambdaFactor(const double value) {}

    virtual double getIncreaseFactor() {}
    virtual double getDecreaseFactor() {}
    virtual double getInitLambdaFactor() {}

protected:
    size_t iteration_times;       // recent iteration times
    size_t max_iteration_times;  // maximum iteration times

    double upper_bound; // upper bound for stepsize
    double lower_bound; // lower bound for stepsize
};
}

#endif
