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
    using AccurateBase::lambda;

    using AccurateBase::f;

public:
    /// @brief Constructors
    DichotomousMethod() { params = new AccurateSearchParams; }
    DichotomousMethod(AccurateSearchParams* given_params) { params = given_params; }

    ~DichotomousMethod() { delete params; }

    void printProcess()
    {
        std::cout << "interative times: " << params->getIterationTimes()
                  << "  " << "[alpha, lambda, mu, beta]: "
                  << "[" << alpha << ", " << lambda << ", "
                  << mu << ", " << beta << "]" << std::endl;
    }

    T search(JacobianType& d) override
    {
        this->reset(params);
        auto epsilon = params->getStepsizeAccuracy();
        lambda = (alpha + beta)/2 - epsilon;
        mu = (alpha + beta)/2 + epsilon;

        params->setMaxIterations(int(std::log2((beta - alpha)/epsilon)));

        while (true)
        {
            //this->printProcess();
            // stoping condition: (alpha - beta) < 2*epsilon or reach max iteration times
            if (params->getIterationTimes() == params->getMaxIterations())
            {
                auto result = lambda;
                this->reset(params);
                return result;
            }
            params->nextIteration();
            // if f(x + lambda * d) < f(x + mu * d)
            if ((*f)(f->getX() + lambda * d.transpose()) < (*f)(f->getX() + mu * d.transpose()))
            {
                // cut the right half of the interval
                beta = mu;
                lambda = (alpha + beta)/2 - epsilon;
                mu = (alpha + beta)/2 + epsilon;
            }
            else
            {
                // cut the left half of the interval
                alpha = lambda;
                lambda = (alpha + beta)/2 - epsilon;
                mu = (alpha + beta)/2 + epsilon;
            }
        }
    }

protected:
    T mu;       // observation variable
    AccurateSearchParams* params;
};

}

#endif
