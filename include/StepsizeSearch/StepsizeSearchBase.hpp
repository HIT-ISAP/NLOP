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

    /// @brief Set max iteration times
    void setMaxIterations(int value)
    {
        max_iteration_times = value;
    }

    virtual void init(FunctorType* f) = 0;

    FunctorType* f; // Target function

    int iteration_times = 0;
    int max_iteration_times = 20;
};
}

#endif
