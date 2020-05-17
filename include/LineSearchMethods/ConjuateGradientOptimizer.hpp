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

    using JacobianType = typename Functor<T, N>::JacobianType;

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

        stepsize_searcher = new GoldSectionMethod<T, PhiFunctortype>;
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
        while (true){
            this->updateValueAndJacobian();
            if (cg_params->iteration_times > cg_params->max_iteration_times)
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                this->printResult();
                return f->getX();
            }
            if (f->getJacobian().norm() < cg_params->min_gradient)
            {
                std::cout << "Iteration times: " << cg_params->iteration_times << std::endl;
                this->printResult();
                return f->getX();
            }
            else
            {
                cg_params->iteration_times++;

                stepsize_searcher->init();

                g = f->getJacobian();
                if (cg_params->iteration_times == 1)
                {
                    beta = 0;
                    last_d.setZero(1, N);
                    last_g.setZero(1, N);
                }
                else
                {
                    beta = (g*g.transpose()/(last_g*last_g.transpose()))(0,0);
                }

                d = -g + beta*last_d;
                stepsize_searcher->phi.setP(f->getX());
                stepsize_searcher->phi.setJ(-d);

                stepsize = stepsize_searcher->search();
                f->setX(f->getX() + stepsize * d.transpose());

                last_g = g;
                last_d = d;
                //this->printProcessInformation();
            }
        }

    }

    ConjuateGradientParams* cg_params;
    JacobianType last_g ;
    JacobianType g;
    JacobianType d;
    JacobianType last_d;

    T beta = 0;
};
}

#endif
