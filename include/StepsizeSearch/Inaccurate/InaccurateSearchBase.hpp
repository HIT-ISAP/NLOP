#ifndef INACCURATESEARCHBASE_HPP
#define INACCURATESEARCHBASE_HPP

#include <Utils/Functor.hpp>
#include <iostream>

namespace NLOP {

/// @class NLOP::InaccurateSearchBase
/// @brief Abstract base class for all inaccurate search methods
/// @param T The numeric scalar type
template <typename FunctorType>
class InaccurateSearchBase
{
public:
    using T = typename FunctorType::Scalar;
    using InputType = typename FunctorType::InputType;
    using JacobianType = typename FunctorType::JacobianType;
    /// @brief Constructor
    InaccurateSearchBase() {}

    /// @brief Iteratively search the inaccurate stepsize
    virtual T search(JacobianType& d) = 0;

    FunctorType* f; // Target function
    T alpha = 1.5; // Increase factor
    T beta = 0.5; // Decrease factor
    T lambda; // Stepsize

    int iteration_times = 0;
    int max_iteration_times = 10;

};
}

#endif
