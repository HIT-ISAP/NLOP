#ifndef GOLDENSECTIONMETHOD_HPP
#define GOLDENSECTIONMETHOD_HPP

#include <StepsizeSearch/Accurate/AccurateSearchBase.hpp>

namespace NLOP {

/// @class NLOP::GoldSectionMethod
/// @brief Golden Section Method (t = 0.618)
/// @param FunctorType Target function
template<typename FunctorType>
class GoldSectionMethod: public AccurateSearchBase<FunctorType>
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
    /// @brief Constructor
    GoldSectionMethod() { params = new AccurateSearchParams; }
    //GoldSectionMethod(StepsizeSearchParamsBase given_params) { *params = given_params; }

    ~GoldSectionMethod() { delete params; }

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
        std::cout << "iteration times: " << params->getIterationTimes()
                  << "  " << "[alpha, lambda, mu, beta]: "
                  << "[" << alpha << ", " << lambda << ", "
                  << mu << ", " << beta << "]" << std::endl;
    }
    /// @brief Use Golden Section Method to search the optimal stepsize
    T search(JacobianType& d) override
    {
        this->reset(params);
        lambda = alpha + (1 - t) * (beta - alpha);
        mu = alpha + t * (beta - alpha);
        while (true)
        {
            //this->printProcess();
            // stoping condition: (alpha - beta) < epsilon
            if ((beta - alpha) < params->getStepsizeAccuracy())
            {
                auto stepsize = (alpha + beta)/2;
                //this->reset(params); // reset alpha, beta for next searching
                return stepsize;
            }
            params->nextIteration();
            // if f(x + lambda * d) > f(x + mu * d)
            if ((*f)(f->getX() + lambda * d.transpose()) > (*f)(f->getX() + mu * d.transpose()))
            {
                // cut the interval [alpha, lambda)
                alpha = lambda;
                lambda = mu;
                mu = alpha + t * (beta - alpha);
            }
            else
            {
                // cut the interval (mu, beta]
                beta = mu;
                mu = lambda;
                lambda = alpha + (1 - t) * (beta - alpha);
            }
        }
    }

private:
    T mu;
    T t = 0.618;                // reduction ratio

    AccurateSearchParams* params;
};

}

#endif
