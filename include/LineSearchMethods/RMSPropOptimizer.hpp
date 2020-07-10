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
    /// @brief Constructor and Deconstructor
    RMSPropOptimizer() {}
    ~RMSPropOptimizer() {}

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
        auto max_iterations = params->getMaxIterations();
        auto min_gradient = params->getMinGradient();
        auto alpha = params->getAlpha();
        auto gamma = params->getGamma();
        auto epsilon = params->getEpsilon();

        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/RMSProp.txt");
        while (true) {
            this->updateValueAndJacobian();
            this->printProcessInformation(params);
            if (this->writer.is_open())
                this->writeInformation();
            if (params->getIterationTimes()> max_iterations)
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
                    s[i] = gamma * s_last[i] + (1 - gamma) * g[i] * g[i];
                    x_next[i] = x[i] - alpha / (epsilon + sqrt(s[i])) * g[i];
                }

                // update x
                f->setX(x_next);
                s_last = s;

                params->nextIteration();
            }
        }
    }

private:
    RMSPropParams* params;
    JacobianType s;         // cumulative sum of squares of gradients
    JacobianType s_last;    // cumulative sum of squares of gradients at last time
    InputType x;            // x at time k
    InputType x_next;       // x at time k+1
    JacobianType g;         // gradient at time k
};
}

#endif
