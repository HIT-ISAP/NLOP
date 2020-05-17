#ifndef STEEPESTDESCENT_HPP
#define STEEPESTDESCENT_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/SteepestDescentParams.hpp>

namespace NLOP {

/// @class SteepestDescentOptimizer
/// @brief Steepest descent method optimizer
/// @param T The numeric scalar type
/// @param N The dimension of variable x
template<typename T, int N, typename FunctorType, typename PhiFunctortype>
class SteepestDescentOptimizer: public LineSearchOptimizer<T, N, FunctorType, PhiFunctortype>
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
    SteepestDescentOptimizer() {}


    ~SteepestDescentOptimizer() {}

    /// @brief Initialize SteepestDescentOptimizer
    void init(const InputType& initial, FunctorType* f, SteepestDescentParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->sd_params = params;

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

    /// @brief Steepest Descent optimization process
    InputType optimize() override
    {
        //sd_params->print();
        std::cout << "Initial Configurations: " << "\n"
                  << "x0: (" << f->getX().transpose() << ") \n"
                  << "f(x0) = " << f->getY() << std::endl;
        // stepsize_searcher->print();
        while (true) {
            this->updateValueAndJacobian();
            if (sd_params->iteration_times > sd_params->max_iteration_times)
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                return f->getX();
            }
            if (f->getJacobian().norm() < sd_params->min_gradient)
            {
                this->printResult();
                return f->getX();
            }
            else
            {
                sd_params->iteration_times++;

                stepsize_searcher->init();

                stepsize_searcher->phi.setP(f->getX());
                stepsize_searcher->phi.setJ(f->getJacobian());

                stepsize = stepsize_searcher->search();
                f->setX(f->getX() + stepsize * (-(f->getJacobian()).transpose()));
                // this->printProcessInformation();
            }
        }
    }

    SteepestDescentParams* sd_params;
};
}

#endif
