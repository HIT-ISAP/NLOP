#ifndef STEPSIZESEARCHBASE_HPP
#define STEPSIZESEARCHBASE_HPP

#include "Utils/Functor.hpp"
#include "Utils/Matrix.hpp"
#include <iostream>

namespace NLOP {
/// @class NLOP::StepsizeSearchBase
/// @brief Abstract base class for all stepsize search methods
/// @param FunctorType Target function
template<typename FunctorType>
class StepsizeSearchBase
{
protected:
    using T = typename FunctorType::Scalar;
    using InputType = typename FunctorType::InputType;
    using ValueType = typename FunctorType::ValueType;
    using JacobianType = typename FunctorType::JacobianType;

public:
    /// @brief Iteratively search the accurate optimal stepsize
    virtual T search(JacobianType& d) = 0;

    /// @brief Set and get max iteration times
    void setMaxIterations(size_t value) { max_iteration_times = value; }
    size_t getMaxIterations() const { return max_iteration_times; }

    /// @brief Set and get recent iteration times
    void setIterationTimes(size_t value) { iteration_times = value; }
    size_t getIterationTimes() const { return iteration_times; }

    virtual void init(FunctorType* f) = 0;

    FunctorType* f; // Target function

    size_t iteration_times = 0;
    size_t max_iteration_times = 20;
};
}

#endif
