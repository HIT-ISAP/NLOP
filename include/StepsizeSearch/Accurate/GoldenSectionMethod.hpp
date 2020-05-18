#ifndef GOLDENSECTIONMETHOD_HPP
#define GOLDENSECTIONMETHOD_HPP

#include <StepsizeSearch/Accurate/AccurateSearchBase.hpp>

namespace NLOP {

/// @class NLOP::GoldSectionMethod
/// @brief Golden Section Method (t = 0.618)
/// @param T The numeric scalar type

template<typename T, typename FunctorType>
class GoldSectionMethod: public AccurateSearchBase<T, FunctorType>
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
    /// @brief Use Golden Section Method to search the optimal stepsize
    T search() override
    {
        double t = 0.618; // reduction ratio
        lambda = alpha + (1 - t) * (beta - alpha);
        mu = alpha + t * (beta - alpha);
        while (true) {

            /*
            std::cout << "interative times: " << iteration_times << "  "
                 << "  " << "[alpha, lambda, mu, beta]: "
                 << "[" << alpha << ", " << lambda << ", "
                 << mu << ", " << beta << "]" << std::endl;
                 */

            /// stoping condition: (alpha - beta) < epsilon
            if ((beta - alpha) < epsilon)
            {
                //std::cout << "Finished! Optimal lambda: "
                //          << (alpha + beta)/2 << std::endl;

                return (alpha + beta)/2;
            }
            iteration_times++;

            if (phi(lambda) > phi(mu))
            {
                /// cut the interval [alpha, lambda)
                alpha = lambda;
                lambda = mu;
                mu = alpha + t * (beta - alpha);
            }
            else
            {
                /// cut the interval (mu, beta]
                beta = mu;
                mu = lambda;
                lambda = alpha + (1 - t) * (beta - alpha);
            }
        }
    }

private:
    T mu;

};

}

#endif
