#ifndef GOLDSTEINMETHOD_HPP
#define GOLDSTEINMETHOD_HPP

#include <StepsizeSearch/Inaccurate/InaccurateSearchBase.hpp>
#include <StepsizeSearchParams/GoldsteinParams.hpp>

namespace NLOP {

template<typename FunctorType>
class GoldsteinMethod: public InaccurateSearchBase<FunctorType>
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
    /// @brief Constructor and Deconstructor
    GoldsteinMethod() { params = new GoldsteinParams; }
    ~GoldsteinMethod() { delete params; }

    /// @brief Set params
    void setParams(StepsizeSearchParamsBase* given_params) override
    {
        params->setUpperBound(given_params->getUpperBound());
        params->setLowerBound(given_params->getLowerBound());
        params->setMaxIterations(given_params->getMaxIterations());
        params->setIncreaseFactor(given_params->getIncreaseFactor());
        params->setDecreaseFactor(given_params->getDecreaseFactor());
    }

    /// @brief Search inaccurate stepsize iteratively using Goldstein method
    /// @param d Direction for stepsize searching
    T search(JacobianType& d) override
    {
        // get params
        auto max_iteration_times = params->getMaxIterations();
        auto alpha = params->getIncreaseFactor();
        auto beta = params->getDecreaseFactor();
        auto rho = params->getRho();

        // initialization
        this->reset(params);
        lambda = params->getInitLambdaFactor() * params->getUpperBound();

        while (true)
        {
            if (params->getIterationTimes() > max_iteration_times)
            {
                std::cerr << "Beyong the max iteration times!" << std::endl;
                return lambda;
            }
            params->nextIteration();
            lhs = (*f)(f->getX() + lambda * d.transpose()) - f->getY();
            rhs1 = (rho * f->getJacobian() * d.transpose() * lambda)(0,0);
            // condition (1): f(x(k+1)) - f(x(k)) <= rho * J(x(k)) * lambda * d
            // condition (2): f(x(k+1)) - f(x(k)) >= (1 - rho) * J(x(k)) * lambda * d
            // if condition (1) is satisfied
            if (lhs <= rhs1)
            {
                rhs2 = ((1-rho) * f->getJacobian() * d.transpose() * lambda)(0,0);
                // if condition (1) and (2) are satisfied simultaneously, stop
                if (lhs >= rhs2)
                    return lambda;
                else
                    lambda *= alpha; // increase Stepsize
            }
            else
                lambda *= beta; // decrease Stepsize
        }
    }

private:
    T lhs;      // left hand side of condition (1) and (2)
    T rhs1;     // right hand side of condition (1)
    T rhs2;     // right hand side of condition (2)

    GoldsteinParams* params;
};
}


#endif
