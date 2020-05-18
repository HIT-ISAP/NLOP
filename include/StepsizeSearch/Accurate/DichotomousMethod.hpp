#ifndef DICHOTOMOUSMETHOD_HPP
#define DICHOTOMOUSMETHOD_HPP

#include <StepsizeSearch/Accurate/AccurateSearchBase.hpp>

namespace NLOP {

/// @class NLOP::Dichotomous Method
/// @brief Dichotomous Method for General Continues Functions
/// @param T the numberic scalar type
template<typename T, typename FunctorType>
class DichotomousMethod: public AccurateSearchBase<T, FunctorType>
{
protected:
    using AccurateSearchBase<T, FunctorType>::alpha;
    using AccurateSearchBase<T, FunctorType>::beta;
    using AccurateSearchBase<T, FunctorType>::epsilon;
    using AccurateSearchBase<T, FunctorType>::lambda;
    using AccurateSearchBase<T, FunctorType>::iteration_times;
    using AccurateSearchBase<T, FunctorType>::max_iteration_times;
    using AccurateSearchBase<T, FunctorType>::phi;
public:
    T search() override
    {
        while (true)
        {
            lambda = (alpha + beta)/2 - epsilon;
            mu = (alpha + beta)/2 + epsilon;
            std::cout << "interative times: " << iteration_times
                      << "  " << "[alpha, lambda, mu, beta]: "
                      << "[" << alpha << ", " << lambda << ", "
                      << mu << ", " << beta << "]" << std::endl;

            /// stoping condition: (alpha - beta) < 2*epsilon or reach max iteration times
            if (beta - alpha < 2*epsilon || iteration_times == max_iteration_times)
            {
                std::cout << "Finished! Optimal lambda: "
                          << lambda << std::endl;
                return lambda;
            }
            iteration_times++;
            if (phi(mu) > phi(lambda))
            {
                /// cut the right half of the interval
                beta = mu;
                lambda = (alpha + beta)/2 - epsilon;
                mu = (alpha + beta)/2 + epsilon;
            }
            else
            {
                /// cut the left half of the interval
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
