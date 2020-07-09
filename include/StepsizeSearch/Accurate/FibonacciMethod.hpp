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
    using AccurateBase::lambda;

    using AccurateBase::f;

public:
    /// @brief Constructors
    FibonacciMethod() { params = new AccurateSearchParams; }
    //FibonacciMethod(StepsizeSearchParamsBase* given_params) { params = given_params; }
    ~FibonacciMethod() { delete params; }

    /// @brief Set params
    void setParams(StepsizeSearchParamsBase* given_params) override
    {
        params->setUpperBound(given_params->getUpperBound());
        params->setLowerBound(given_params->getLowerBound());
        params->setStepsizeAccuracy(given_params->getStepsizeAccuracy());
        params->setMaxIterations(given_params->getMaxIterations());
    }

    /// @brief Print stepsize searching process information
    void printProcess()
    {
        std::cout << "interative times: " << params->getIterationTimes()
                  << "  " << "[alpha, lambda, mu, beta]: "
                  << "[" << alpha << ", " << lambda << ", "
                  << mu << ", " << beta << "]" << std::endl;
    }

    /// @brief search stepsize using
    T search(JacobianType& d) override
    {
        this->reset(params);
        n = params->getMaxIterations();
        t = fibonacci(n-1)/fibonacci(n);
        lambda = alpha + (1 - t) * (beta - alpha);
        mu = alpha + t * (beta - alpha);
        while (true)
        {
            //this->printProcess();
            if (params->getIterationTimes() > params->getMaxIterations())
            {
                std::cout << "Beyong max iteration times" << std::endl;
                this->reset(params);
                return (beta + alpha)/2;
            }
            // stopping condition: (alpha - beta) < epsilon
            if (beta - alpha < params->getStepsizeAccuracy())
            {
                auto result = (beta + alpha)/2;
                this->reset(params);
                return result;
            }
            params->nextIteration();
            // if f(x + lambda * d) > f(x + mu * d)
            if ((*f)(f->getX() + lambda * d.transpose()) > (*f)(f->getX() + mu * d.transpose()))
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
            // update reduction ratio;
            t = fibonacci(n-1)/fibonacci(n);
        }
    }

private:
    size_t n;   // the expected index of Fibonacci array (default = 13)
    T mu;       // observation variable
    T t;        // reduction ratio
    AccurateSearchParams* params;

};
}

#endif
