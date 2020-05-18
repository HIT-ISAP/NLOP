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
    /// @brief Use Golden Section Method to search the optimal stepsize
    T search(JacobianType& d) override
    {
        //double t = 0.618; // reduction ratio
        lambda = alpha + (1 - t) * (beta - alpha);
        mu = alpha + t * (beta - alpha);
        while (true)
        {
            //this->printProcess();
            // Stoping condition: (alpha - beta) < epsilon
            if ((beta - alpha) < epsilon)
            {
                auto stepsize = (alpha + beta)/2;
                //std::cout << "Finished! Optimal lambda: " << stepsize << std::endl;
                this->reset(); // Reset alpha, beta for next searching
                return stepsize;
            }
            iteration_times++;

            // If f(x + lambda*d) > f(x + mu*d)
            if ((*f)(f->getX()+lambda*d.transpose()) > (*f)(f->getX()+mu*d.transpose()))
            {
                // Cut the interval [alpha, lambda)
                alpha = lambda;
                lambda = mu;
                mu = alpha + t * (beta - alpha);
            }
            else
            {
                // Cut the interval (mu, beta]
                beta = mu;
                mu = lambda;
                lambda = alpha + (1 - t) * (beta - alpha);
            }
        }
    }

private:
    T mu;
    T t = 0.618; // reduction ratio
};

}

#endif
