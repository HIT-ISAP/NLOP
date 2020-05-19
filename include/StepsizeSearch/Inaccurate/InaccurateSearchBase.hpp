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

    using Base::iteration_times;
    using Base::max_iteration_times;
    using Base::f;

    void reset()
    {
        iteration_times = 0;
    }

    void printResult()
    {
        std::cout << "Inaccurate stepsize: " << lambda << std::endl;
        std::cout << "Iteration times: " << iteration_times << std::endl;
    }

    //FunctorType* f; // Target function
    T alpha = 1.5; // Increase factor
    T beta = 0.5; // Decrease factor
    T lambda; // Stepsize
};
}

#endif
