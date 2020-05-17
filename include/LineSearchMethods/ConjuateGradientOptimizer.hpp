#ifndef CONJUATEGRADIENT_HPP
#define CONJUATEGRADIENT_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/ConjuateGradientParams.hpp>

namespace NLOP {
/// @class ConjuateGradientOptimizer
/// @brief Conjuate gradient method optimizer
/// @param T The numeric scalar type
/// @param N The dimension of variable x
template<typename T, int N, typename FunctorType, typename PhiFunctortype>
class ConjuateGradientOptimizer: public LineSearchOptimizer<T, N, FunctorType, PhiFunctortype>
{
protected:
    using LineSearch = LineSearchOptimizer<T, N, FunctorType, PhiFunctortype>;

    using LineSearch::f;
    using LineSearch::params;
    using LineSearch::stepsize_searcher;
    using LineSearch::stepsize;
    using typename LineSearch::InputType;

public:
    /// @brief Constructors
    ConjuateGradientOptimizer() {}

    ~ConjuateGradientOptimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f, ConjuateGradientParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->cg_params = params;

        // stepsize_searcher = new GoldSectionMethod<T, PhiFunctortype>;
        /*
        /// Initial stepsize search method
        switch (this->params->stepsize_search_method) {
        case SteepestDescentParams::GOLDENSECTION:
            stepsize_searcher = new GoldSectionMethod<T, FunctorType>;
            break;
        case SteepestDescentParams::DICHOTOMOUS:
            stepsize_searcher = new DichotomousMethod<T, FunctorType>;
            break;
        case SteepestDescentParams::NEWTON:
            /// TODO
            break;
        case SteepestDescentParams::FIBONACCI:
            stepsize_searcher = new FibonacciMethod<T, FunctorType>;
            break;
        default:
            break;
        }
        */
    }

    /// @brief Conjuate Gradient optimization process
    InputType optimize() override
    {
        std::cout << "Initial Configurations: " << "\n"
                  << "x0: (" << f->getX().transpose() << ") \n"
                  << "f(x0) = " << f->getY() << std::endl;

    }

    ConjuateGradientParams* cg_params;

};
}

#endif
