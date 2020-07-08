#ifndef WOLFEPOWELLMETHOD_HPP
#define WOLFEPOWELLMETHOD_HPP

#include <StepsizeSearch/Inaccurate/InaccurateSearchBase.hpp>
#include <StepsizeSearchParams/WolfePowellParams.hpp>

namespace NLOP {

template<typename FunctorType>
class WolfePowellMethod: public InaccurateSearchBase<FunctorType>
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
    WolfePowellMethod() { params = new WolfePowellParams; }
    WolfePowellMethod(WolfePowellParams* given_params) { params = given_params; }

    ~WolfePowellMethod() { delete params; }

    /// @brief Compute y and jacobian with given x
    void computeValueAndJacobian(const InputType& x, ValueType* v, JacobianType* jac)
    {
        adjac(x, v, jac);
    }

    /// @brief Search inaccurate stepsize iteratively using Goldstein method
    /// @param d Direction for stepsize searching
    T search(JacobianType& d) override
    {
        this->reset(params);
        // get params
        auto max_iteration_times = params->getMaxIterations();
        auto alpha = params->getIncreaseFactor();
        auto beta = params->getDecreaseFactor();
        auto rho = params->getRho();
        auto sigma = params->getSigma();
        // initial stepsize
        lambda = params->getInitLambdaFactor() * params->getUpperBound();
        while (true)
        {
            if (params->getIterationTimes() > max_iteration_times)
            {
                std::cerr << "Beyong the max iteration times!" << std::endl;
                return lambda;
            }
            params->nextIteration();
            x_next = f->getX() + lambda * d.transpose(); // x(k+1)
            this->computeValueAndJacobian(x_next, &y_next, &jac_next);

            lhs1 = (*f)(x_next) - f->getY();
            rhs1 = (rho * f->getJacobian() * d.transpose() * lambda)(0,0);
            // condition (1): f(x(k+1)) - f(x(k)) <= rho * J(x(k)) * lambda * d
            // condition (2): J(x(k+1)) * lambda * d >= sigma * J(x(k)) * lambda * d
            // if condition (1) is satisfied
            if (lhs1 <= rhs1)
            {
                lhs2 = (jac_next * d.transpose() * lambda)(0,0);
                rhs2 = (sigma * f->getJacobian() * d.transpose() * lambda)(0,0);
                // if condition (1) and (2) are satisfied simultaneously, stop
                if (lhs2 >= rhs2)
                {
                    this->reset(params);
                    return lambda;
                }
                else
                {
                    // Increase Stepsize
                    lambda *= alpha;
                }
            }
            else
            {
                // Decrease Stepsize
                lambda *= beta;
            }
        }
    }

private:
    WolfePowellParams* params;

    JacobianType jac_next;
    ValueType y_next;
    InputType x_next;

    Eigen::AutoDiffJacobian<FunctorType> adjac;

    T lhs1;     // left hand side of condition (1)
    T lhs2;     // left hand side of condition (2)

    T rhs1;     // right hand side of condition (1)
    T rhs2;     // right hand side of condition (2)
};
}

#endif
