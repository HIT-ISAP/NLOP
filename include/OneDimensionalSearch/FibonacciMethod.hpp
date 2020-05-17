#ifndef FIBONACCIMETHOD_HPP
#define FIBONACCIMETHOD_HPP

#include <OneDimensionalSearch/OneDimensionalSearchMethods.hpp>
#include <Utils/Utils.hpp>

namespace NLOP {

/// @class NLOP::FibonacciMethod
/// @brief Fibonacci method (t = 0.618)
/// @param T The numeric scalar type
template<typename T, typename FunctorType>
class FibonacciMethod: public OneDimSearch<T, FunctorType>
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
        double t = fibonacci(n-1)/fibonacci(n);
        lambda = alpha + (1 - t) * (beta - alpha);
        mu = alpha + t * (beta - alpha);
        while (true)
        {
            std::cout << "interative times: " << iteration_times
                      << "  " << "[alpha, lambda, mu, beta]: "
                      << "[" << alpha << ", " << lambda << ", "
                      << mu << ", " << beta << "]" << std::endl;

            /// stoping condition: (alpha - beta) < epsilon
            if (beta - alpha < epsilon)
            {
                std::cout << "Finished! Optimal lambda: "
                          << (beta + alpha)/2 << std::endl;
                return (beta + alpha)/2;
            }
            iteration_times++;

            if (phi(lambda) > phi(mu))
            {
                alpha = lambda;
                lambda = mu;
                mu = alpha + t * (beta - alpha);
            }
            else
            {
                beta = mu;
                mu = lambda;
                lambda = alpha + (1 - t) * (beta - alpha);
            }
            if (--n == 0)
            {
                std::cerr << "n is too small, please choose a bigger n and try again" << std::endl;
                return -1;
            }
            /// update reduction ratio;
            t = fibonacci(n-1)/fibonacci(n);
        }
    }

private:
    /// The expected index of Fibonacci array (default = 13)
    int n = 15;

    T mu;

};
}

#endif
