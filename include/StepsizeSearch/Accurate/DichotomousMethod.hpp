#ifndef DICHOTOMOUSMETHOD_HPP
#define DICHOTOMOUSMETHOD_HPP

#include <StepsizeSearch/Accurate/AccurateSearchBase.hpp>

namespace NLOP {

/// @class NLOP::Dichotomous Method
/// @brief Dichotomous Method for General Continues Functions
/// @param FunctorType Target function
template<typename FunctorType>
class DichotomousMethod: public AccurateSearchBase<FunctorType>
{
protected:
    using AccurateBase = AccurateSearchBase<FunctorType>;
    using typename AccurateBase::T;
    using typename AccurateBase::InputType;
    using typename AccurateBase::ValueType;
    using typename AccurateBase::JacobianType;

    using AccurateBase::alpha;
    using AccurateBase::beta;
    using AccurateBase::epsilon;
    using AccurateBase::lambda;

    using AccurateBase::iteration_times;
    using AccurateBase::max_iteration_times;
    using AccurateBase::f;

public:
    void printProcess() override
    {
        std::cout << "interative times: " << iteration_times
                  << "  " << "[alpha, lambda, mu, beta]: "
                  << "[" << alpha << ", " << lambda << ", "
                  << mu << ", " << beta << "]" << std::endl;
    }

    T search(JacobianType& d) override
    {
        epsilon = 0.0005;
        lambda = (alpha + beta)/2 - epsilon;
        mu = (alpha + beta)/2 + epsilon;

        max_iteration_times = int(std::log2((beta - alpha)/epsilon));

        while (true)
        {
            //this->printProcess();

            // Stoping condition: (alpha - beta) < 2*epsilon or reach max iteration times
            if (iteration_times == max_iteration_times)
            {
                auto result = lambda;
                //std::cout << "Finished! Optimal lambda: " << result << std::endl;
                this->reset();
                return result;
            }
            iteration_times++;

            // If f(x + lambda*d) < f(x + mu*d)
            if ((*f)(f->getX()+lambda*d.transpose()) < (*f)(f->getX()+mu*d.transpose()))
            {
                // Cut the right half of the interval
                beta = mu;
                lambda = (alpha + beta)/2 - epsilon;
                mu = (alpha + beta)/2 + epsilon;
            }
            else
            {
                // Cut the left half of the interval
                alpha = lambda;
                lambda = (alpha + beta)/2 - epsilon;
                mu = (alpha + beta)/2 + epsilon;
            }
        }
    }

protected:
    T mu;
};

}

#endif
