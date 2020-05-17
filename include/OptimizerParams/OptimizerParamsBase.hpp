#ifndef OPTIMIZERPARAMSBASE_HPP
#define OPTIMIZERPARAMSBASE_HPP

#include <Utils/Functor.hpp>

namespace NLOP {

/// @class NLOP::OptimizerParamsBase
/// @brief Abstract base class for params of all non-linear optimization methods

class OptimizerParamsBase
{
public:
    /// @brief Use default optimizer params
    // virtual void setDefaults() = 0;

    /// @brief Set max iteration times
    void setMaxItertaions(int max)
    {
        max_iteration_times = max;
    }

    /// @brief set verbosity
    void setVerbosity(bool value)
    {
        verbosity = value;
    }

    /// @brief Print optimizer params
    //virtual void print();

    //virtual ~OptimizerParamsBase() {}

    int iteration_times = 0;
    int max_iteration_times = 10000; // Default = 10000;
    bool verbosity = false; // Default = false;

};
}


#endif
