#ifndef CONSTRAINED_OPTIMIZER_BASE_HPP
#define CONSTRAINED_OPTIMIZER_BASE_HPP

#include <OptimizerBase/OptimizerBase.hpp>
#include <OptimizerParams/OptimizerParamsBase.hpp>
#include <StepsizeSearch/StepsizeSearchBase.hpp>
#include <Utils/Constraint.hpp>

namespace NLOP {
/// @class NLOP::ConstrainedOptimizerBase
/// @brief Abstract base class for all constrained optimization methods
/// @param FunctorType Target function type
/// @param ConstraintType Constraint type
template<typename FunctorType, typename ConstraintType>
class ConstrainedOptimizerBase
{
protected:
    using T = typename FunctorType::Scalar;
    using InputType = typename FunctorType::InputType;
    using ValueType = typename FunctorType::ValueType;
    using ResidualType = typename ConstraintType::ValueType;
    using JacobianType = typename FunctorType::JacobianType;
    using HessianType = typename FunctorType::HessianType;

public:
    struct AugmentedLagrangianFunctor: public NLOP::Functor<T, f->m_inputs>
    {
        void operator() (const InputType& x, ValueType* v, JacobianType* _j = 0) const
        {
            (*v)[0] = f->getY() + lambda.dot(e) + rho * e.dot(e);
        }

        void operator ()(const ActiveInput& x, ActiveValue* v) const
        {
            (*v)[0] = (*f)(x) + lambda.dot((h->A)*x - h->b) + rho * ((h->A)*x - h->b).dot(((h->A)*x - h->b));
        }
    };

protected:
     ConstraintType* h;
     FunctorType* f;
     AugmentedLagrangianFunctor* l;

     ResidualType lambda;
     ResidualType e;
     T rho;

     // Tool to compute value and jacobian
     Eigen::AutoDiffJacobian<AugmentedLagrangianFunctor> adjac;
};
}

#endif
