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
        this->printInitialConfigurations();
        this->writer.open("../data/"
                          "Adagrad.txt");
        while (true){
            this->updateValueAndJacobian();
            this->writeInformation();
            if (params->iteration_times > params->max_iteration_times)
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                this->printResult();
                return f->getX();
            }
            if (f->getJacobian().norm() < params->min_gradient)
            {
                std::cout << "Iteration times: " << params->iteration_times << std::endl;
                this->printResult();
                return f->getX();
            }
            else
            {
                params->iteration_times++;
                this->printProcessInformation();

                x = f->getX();
                g = f->getJacobian();

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                {
                    s[i] += g[i] * g[i];
                    x_next[i] = x[i] - params->alpha / (params->epsilon + sqrt(s[i])) * g[i];
                }

                f->setX(x_next);
            }
        }
        this->writer.close();
    }


private:
    AdagradParams* params;
    JacobianType s; // cumulative sum of squares of gradients
    JacobianType g; // gradient at time k
    InputType x_next; // x at time k+1
    InputType x; // x at time k
};
}

#endif
