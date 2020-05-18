#ifndef GOLDSTEINMETHOD_HPP
#define GOLDSTEINMETHOD_HPP

#include <StepsizeSearch/Inaccurate/InaccurateSearchBase.hpp>

namespace NLOP {

template<typename FunctorType>
class GoldsteinMethod: public InaccurateSearchBase<FunctorType>
{
protected:
    using InaccurateSearchBase<FunctorType>::alpha;
    using InaccurateSearchBase<FunctorType>::beta;
    using InaccurateSearchBase<FunctorType>::lambda;
    using InaccurateSearchBase<FunctorType>::iteration_times;
    using InaccurateSearchBase<FunctorType>::max_iteration_times;
    using InaccurateSearchBase<FunctorType>::f;
    using typename InaccurateSearchBase<FunctorType>::T;
    using typename InaccurateSearchBase<FunctorType>::InputType;
    using typename InaccurateSearchBase<FunctorType>::JacobianType;

public:

    /// @brief Set rho
    void setRho(T new_value)
    {
        rho = new_value;
    }

    /// @brief Initialization for stepsize search
    /// @param f Target function
    /// @param alpha Increase factor
    /// @param beta Decrease factor
    void init(FunctorType* f, T alpha = 1.5, T beta = 0.5, T max_lambda = 1)
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

            lhs = (*f)(f->getX()+lambda*d.transpose()) - f->getY();

            rhs1 = (rho * f->getJacobian() * d.transpose() * lambda)(0,0);

            // If condition (1) is satisfied
            if (lhs <= rhs1)
            {
                rhs2 = ((1-rho) * f->getJacobian() * d.transpose() * lambda)(0,0);

                // If condition (1) and (2) are satisfied simultaneously, stop
                if (lhs >= rhs2)
                {
                    // std::cout << "Inaccurate stepsize: " << lambda << std::endl;
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
    T rho = 0.1;
    T max_lambda = 1; // max stepsize

    T init_lambda_factor = 0.1; // initial stepsize factor

    T lhs; // left hand side of condition (1) and (2)
    T rhs1; // right hand side of condition (1)
    T rhs2; // right hand side of condition (2)
};
}


#endif
