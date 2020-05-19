#ifndef ARMIJOMETHOD_HPP
#define ARMIJOMETHOD_HPP

#include <StepsizeSearch/Inaccurate/InaccurateSearchBase.hpp>

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

    using InaccurateBase::alpha;
    using InaccurateBase::beta;
    using InaccurateBase::lambda;

    using InaccurateBase::iteration_times;
    using InaccurateBase::max_iteration_times;
    using InaccurateBase::f;

public:
    /// @brief Constructor
    ArmijoMethod() {}

    /// @brief Set rho
    void setRho(T new_value)
    {
        rho = new_value;
    }

    /// @brief Set mu
    void setMu(T new_value)
    {
        mu = new_value;
    }

    /// @brief Initialization for stepsize search
    /// @param f Target function
    void init(FunctorType* f) override
    {
        this->f = f;
        this->alpha = alpha;
        this->beta = beta;
        this->max_lambda = max_lambda;
    }

    /// @brief Search inaccurate stepsize iteratively using Goldstein method
    /// @param d Direction for stepsize searching
    T search(JacobianType& d) override
    {
        // Initial stepsize
        lambda = init_lambda_factor * max_lambda;
        while (true)
        {
            if (iteration_times > max_iteration_times)
            {
                std::cerr << "Beyong the max iteration times!" << std::endl;
                return lambda;
            }

            iteration_times++;

            lhs = (*f)(f->getX()+lambda*d.transpose()) - f->getY();

            rhs1 = (rho * f->getJacobian() * d.transpose() * lambda)(0,0);

            // Condition (1): f(x(k+1)) - f(x(k)) <= rho * J(x(k)) * lambda * d
            // Condition (2): f(x(k+1)) - f(x(k)) >= (mu * rho) * J(x(k)) * lambda * d

            // If condition (1) is satisfied
            if (lhs <= rhs1)
            {
                rhs2 = ((mu * rho) * f->getJacobian() * d.transpose() * lambda)(0,0);

                // If condition (1) and (2) are satisfied simultaneously, stop
                if (lhs >= rhs2)
                {
                    //this->printResult();
                    this->reset();
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
    T rho = 0.1; // rho <= 0.1
    T mu = 10; // mu = 5 ~ 10

    T max_lambda = 1; // max stepsize
    T init_lambda_factor = 0.1; // initial stepsize factor

    T lhs; // left hand side of condition (1) and (2)
    T rhs1; // right hand side of condition (1)
    T rhs2; // right hand side of condition (2)
};
}

#endif
