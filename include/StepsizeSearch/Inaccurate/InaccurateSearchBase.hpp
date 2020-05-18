#ifndef INACCURATESEARCHBASE_HPP
#define INACCURATESEARCHBASE_HPP

#include <StepsizeSearch/StepsizeSearchBase.hpp>

namespace NLOP {

/// @class NLOP::InaccurateSearchBase
/// @brief Abstract base class for all inaccurate search methods
/// @param FunctorType Target function type
template <typename FunctorType>
class InaccurateSearchBase: public StepsizeSearchBase<FunctorType>
{
public:
    using Base = StepsizeSearchBase<FunctorType>;
    using typename Base::T;
    using typename Base::InputType;
    using typename Base::ValueType;
    using typename Base::JacobianType;

    using typename Base::iteration_times;
    using typename Base::max_iteration_times;
    using typename Base::f;

    /// @brief Constructor
    // InaccurateSearchBase() {}

    /// @brief Iteratively search the inaccurate stepsize
    // virtual T search(JacobianType& d) = 0;

    //FunctorType* f; // Target function
    T alpha = 1.5; // Increase factor
    T beta = 0.5; // Decrease factor
    T lambda; // Stepsize

    //int iteration_times = 0;
    //int max_iteration_times = 10;

};
}

#endif
