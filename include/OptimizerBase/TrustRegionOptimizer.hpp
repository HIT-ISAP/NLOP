#ifndef TRUSTREGIONOPTIMIZER_HPP
#define TRUSTREGIONOPTIMIZER_HPP

#include <OptimizerBase/OptimizerBase.hpp>

namespace NLOP {
/// @class NLOP::TrustRegionOptimizer
/// @brief Abstract base class for all trust region methods
/// @param FunctorType Target function type
template<typename FunctorType>
class TrustRegionOptimizer: public OptimizerBase<FunctorType>
{
protected:
    using Base = OptimizerBase<FunctorType>;

    using typename Base::T;

    using typename Base::InputType;
    using typename Base::ValueType;
    using typename Base::JacobianType;
    using typename Base::HessianType;

    using Base::f;

public:
    virtual ~TrustRegionOptimizer() {}

};
}

#endif
