#ifndef DICHOTOMOUSMETHOD_HPP
#define DICHOTOMOUSMETHOD_HPP

#include <OneDimensionalSearch/OneDimensionalSearchMethods.hpp>

namespace NLOP {

/// @class NLOP::Dichotomous Method
/// @brief Dichotomous Method for General Continues Functions
/// @param T the numberic scalar type
template<typename T, typename FunctorType>
class DichotomousMethod: public OneDimSearch<T, FunctorType>
{
protected:
    using OneDimSearch<T, FunctorType>::alpha;
    using OneDimSearch<T, FunctorType>::beta;
    using OneDimSearch<T, FunctorType>::epsilon;
    using OneDimSearch<T, FunctorType>::lambda;
    using OneDimSearch<T, FunctorType>::iteration_times;
    using OneDimSearch<T, FunctorType>::max_iteration_times;
    using OneDimSearch<T, FunctorType>::phi;
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
