#ifndef RMSPROPOPTIMIZER_HPP
#define RMSPROPOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/RMSPropParams.hpp>

namespace NLOP {
/// @class RMSPropOptimizer
/// @brief RMSProp method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class RMSPropOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;

public:
    RMSPropOptimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              RMSPropParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;

        s_last.setZero(1, InputType::RowsAtCompileTime);

    }

    /// @brief RMSProp optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations();
        while (true) {
            this->updateValueAndJacobian();
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

                //s = params->gamma * s_last + (1-params->gamma) * (g.dot(g));

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                {
                    s[i] = params->gamma * s_last[i] + (1-params->gamma) * g[i] * g[i];
                    //s[i] += g[i] * g[i];
                    x_next[i] = x[i] - params->alpha / (params->epsilon + sqrt(s[i])) * g[i];
                }

                f->setX(x_next);
                s_last = s;
            }
        }
    }

private:
    RMSPropParams* params;
    JacobianType s; // cumulative sum of squares of gradients
    JacobianType s_last; // cumulative sum of squares of gradients at last time
    InputType x; // x at time k
    InputType x_next; // x at time k+1
    JacobianType g; // gradient at time k
};
}

#endif
