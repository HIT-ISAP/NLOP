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

    using Base::f;

    void reset(StepsizeSearchParamsBase* params) override { params->setIterationTimes(0); }

protected:
    T lambda;       // stepsize
};
}

#endif
