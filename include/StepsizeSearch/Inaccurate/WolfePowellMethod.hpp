#ifndef WOLFEPOWELLMETHOD_HPP
#define WOLFEPOWELLMETHOD_HPP

#include <StepsizeSearch/Inaccurate/InaccurateSearchBase.hpp>

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

    using InaccurateBase::alpha;
    using InaccurateBase::beta;
    using InaccurateBase::lambda;

    using InaccurateBase::iteration_times;
    using InaccurateBase::max_iteration_times;
    using InaccurateBase::f;

public:
    WolfePowellMethod() {}

    /// @brief Initialization for stepsize search
    /// @param f Target function
    void init(FunctorType* f) override
    {
        this->f = f;
        this->alpha = alpha;
        this->beta = beta;
        this->max_lambda = max_lambda;
    }

    /// @brief Compute y and jacobian with given x
    void computeValueAndJacobian(const InputType& x, ValueType* v, JacobianType* jac)
    {
        adjac(x, v, jac);
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

            x_next = f->getX() + lambda * d.transpose(); // x(k+1)
            this->computeValueAndJacobian(x_next, &y_next, &jac_next);

            lhs1 = (*f)(x_next) - f->getY();
            rhs1 = (rho * f->getJacobian() * d.transpose() * lambda)(0,0);

            // Condition (1): f(x(k+1)) - f(x(k)) <= rho * J(x(k)) * lambda * d
            // Condition (2): J(x(k+1)) * lambda * d >= sigma * J(x(k)) * lambda * d

            // If condition (1) is satisfied
            if (lhs1 <= rhs1)
            {
                lhs2 = (jac_next * d.transpose() * lambda)(0,0);
                rhs2 = (sigma * f->getJacobian() * d.transpose() * lambda)(0,0);

                // If condition (1) and (2) are satisfied simultaneously, stop
                if (lhs2 >= rhs2)
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
    T sigma = 0.8; // sigma ~ (rho, 1)

    T max_lambda = 1; // max stepsize
    T init_lambda_factor = 0.1; // initial stepsize factor

    JacobianType jac_next;
    ValueType y_next;
    InputType x_next;

    Eigen::AutoDiffJacobian<FunctorType> adjac;

    T lhs1; // left hand side of condition (1)
    T lhs2; // left hand side of condition (2)

    T rhs1; // right hand side of condition (1)
    T rhs2; // right hand side of condition (2)
};
}

#endif
