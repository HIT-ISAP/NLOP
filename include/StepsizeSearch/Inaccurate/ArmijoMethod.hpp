#ifndef ARMIJOMETHOD_HPP
#define ARMIJOMETHOD_HPP

#include <StepsizeSearch/Inaccurate/InaccurateSearchBase.hpp>
#include <StepsizeSearchParams/ArmijoParams.hpp>

namespace NLOP {

template<typename FunctorType>
class ArmijoMethod: public InaccurateSearchBase<FunctorType>
{
protected:
    using InaccurateBase = InaccurateSearchBase<FunctorType>;
    using typename InaccurateBase::T;
    using typename InaccurateBase::InputType;
    using typename InaccurateBase::ValueType;
    using typename InaccurateBase::JacobianType;

    using InaccurateBase::lambda;
    using InaccurateBase::f;

public:
    /// @brief Constructor
    ArmijoMethod() { params = new ArmijoParams; }
    ArmijoMethod (ArmijoParams* given_params) { params = given_params; }

    ~ArmijoMethod() { delete params; }

    /// @brief Search inaccurate stepsize iteratively using Goldstein method
    /// @param d Direction for stepsize searching
    T search(JacobianType& d) override
    {
        this->reset(params);
        // get params
        auto max_iterations = params->getMaxIterations();
        auto alpha = params->getIncreaseFactor();
        auto beta = params->getDecreaseFactor();
        auto rho = params->getRho();
        auto mu = params->getMu();
        // initial stepsize
        lambda = params->getInitLambdaFactor() * params->getUpperBound();
        while (true)
        {
            if (params->getIterationTimes() > max_iterations)
            {
                std::cerr << "Beyong the max iteration times!" << std::endl;
                return lambda;
            }
            params->nextIteration();
            lhs = (*f)(f->getX() + lambda * d.transpose()) - f->getY();
            rhs1 = (rho * f->getJacobian() * d.transpose() * lambda)(0,0);
            // condition (1): f(x(k+1)) - f(x(k)) <= rho * J(x(k)) * lambda * d
            // condition (2): f(x(k+1)) - f(x(k)) >= (mu * rho) * J(x(k)) * lambda * d
            // if condition (1) is satisfied
            if (lhs <= rhs1)
            {
                rhs2 = ((mu * rho) * f->getJacobian() * d.transpose() * lambda)(0,0);
                // if condition (1) and (2) are satisfied simultaneously, stop
                if (lhs >= rhs2)
                {
                    this->reset(params);
                    return lambda;
                }
                else
                {
                    // increase Stepsize
                    lambda *= alpha;
                }
            }
            else
            {
                // decrease Stepsize
                lambda *= beta;
            }
        }
    }

private:
    T lhs;      // left hand side of condition (1) and (2)
    T rhs1;     // right hand side of condition (1)
    T rhs2;     // right hand side of condition (2)

    ArmijoParams* params;
};
}

#endif
