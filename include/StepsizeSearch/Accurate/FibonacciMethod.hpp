#ifndef FIBONACCIMETHOD_HPP
#define FIBONACCIMETHOD_HPP

#include <StepsizeSearch/Accurate/AccurateSearchBase.hpp>
#include <Utils/Utils.hpp>

namespace NLOP {

/// @class NLOP::FibonacciMethod
/// @brief Fibonacci method (t = 0.618)
/// @param FunctorType Target function type
template<typename FunctorType>
class FibonacciMethod: public AccurateSearchBase<FunctorType>
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

    void reset() override
    {
        alpha = 0;
        beta = 1;
        iteration_times = 0;
        n = 20;
    }

    T search(JacobianType& d) override
    {
        t = fibonacci(n-1)/fibonacci(n);
        lambda = alpha + (1 - t) * (beta - alpha);
        mu = alpha + t * (beta - alpha);
        while (true)
        {
            //this->printProcess();
            if (iteration_times > n)
            {
                std::cout << "Beyong max iteration times" << std::endl;
                this->reset();
                return (beta + alpha)/2;
            }

            // Stopping condition: (alpha - beta) < epsilon
            if (beta - alpha < epsilon)
            {
                auto result = (beta + alpha)/2;
                //std::cout << "Finished! Optimal lambda: " << result << std::endl;
                this->reset();
                return result;
            }
            iteration_times++;

            // If f(x + lambda*d) > f(x + mu*d)
            if ((*f)(f->getX()+lambda*d.transpose()) > (*f)(f->getX()+mu*d.transpose()))
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
                return 0;
            }
            // Update reduction ratio;
            t = fibonacci(n-1)/fibonacci(n);
        }
    }

private:
    int n = 20; // The expected index of Fibonacci array (default = 13)

    T mu; // Observation variable
    T t; // Reduction ratio

};
}

#endif
