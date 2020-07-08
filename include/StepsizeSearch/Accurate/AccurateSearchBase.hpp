#ifndef ACCURATESEARCHBASE_HPP
#define ACCURATESEARCHBASE_HPP

#include <StepsizeSearch/StepsizeSearchBase.hpp>
#include <StepsizeSearchParams/AccurateSearchParams.hpp>

namespace NLOP{
/// @class NLOP::AccurateSearchBase
/// @brief Abstract base class for all accurate search methods
/// @param FunctorType Target function
template<typename FunctorType>
class AccurateSearchBase: public StepsizeSearchBase<FunctorType>
{
protected:
    using Base = StepsizeSearchBase<FunctorType>;
    using typename Base::T;
    using typename Base::InputType;
    using typename Base::ValueType;
    using typename Base::JacobianType;

    using Base::f;

public:
    /*
    /// @brief Initial with default params
    void init(StepsizeSearchParamsBase* params) override
    {
        alpha = params->getLowerBound();
        beta = params->getUpperBound();
        params->setIterationTimes(0);
    }
    */

    /// @brief Reset observation variables and iteration times for the next stepsize searching
    void reset(StepsizeSearchParamsBase* params) override
    {
        alpha = params->getLowerBound();
        beta = params->getUpperBound();
        params->setIterationTimes(0);
    }

protected:
    T alpha;    // observation variables
    T beta;
    T lambda;
};

}

#endif
