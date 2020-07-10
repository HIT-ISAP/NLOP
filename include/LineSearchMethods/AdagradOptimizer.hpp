#ifndef ADAGRADOPTIMIZER_HPP
#define ADAGRADOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/AdagradParams.hpp>

namespace NLOP {
/// @class AdagradOptimizer
/// @brief Adagrad method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class AdagradOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;

public:
    AdagradOptimizer() {}
    ~AdagradOptimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              AdagradParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
    }

    /// @brief Adagrad optimization process
    InputType optimize() override
    {
        // get params
        auto max_iterations = params->getMaxIterations();
        auto min_gradient = params->getMinGradient();
        auto alpha = params->getAlpha();
        auto epsilon = params->getEpsilon();

        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/Adagrad.txt");
        while (true){
            this->updateValueAndJacobian();
            this->printProcessInformation(params);
            if (this->writer.is_open())
                this->writeInformation();
            if (params->getIterationTimes() > max_iterations)
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                if (this->writer.is_open())
                    this->writer.close();
                return f->getX();
            }
            if (f->getJacobian().norm() < min_gradient)
            {
                this->printResult(params);
                if (this->writer.is_open())
                    this->writer.close();
                return f->getX();
            }
            else
            {
                x = f->getX();
                g = f->getJacobian();

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                {
                    s[i] += g[i] * g[i];
                    x_next[i] = x[i] - alpha / (epsilon + sqrt(s[i])) * g[i];
                }

                // update x
                f->setX(x_next);
                params->nextIteration();
            }
        }
    }

private:
    AdagradParams* params;
    JacobianType s;         // cumulative sum of squares of gradients
    JacobianType g;         // gradient at time k
    InputType x_next;       // x at time k+1
    InputType x;            // x at time k
};
}

#endif
