#ifndef ADADELTAOPTIMIZER_HPP
#define ADADELTAOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/AdaDeltaParams.hpp>


namespace NLOP {
/// @class AdaDeltaOptimizer
/// @brief AdaDelta method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class AdaDeltaOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;
    using LineSearch::stepsize;
    using LineSearch::ss;
    using LineSearch::d;

public:
    AdaDeltaOptimizer() {}
    ~AdaDeltaOptimizer() { delete ss; }

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              AdaDeltaParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        this->initStepsizeMethod(params);
        this->ss->bind(this->f);
        this->ss->setParams(params->stepsize_params);

        s.setZero(1, InputType::RowsAtCompileTime);
        s_last.setZero(1, InputType::RowsAtCompileTime);
        s_dx.setZero(InputType::RowsAtCompileTime, 1);
        s_dx_last.setZero(InputType::RowsAtCompileTime, 1);
    }

    /// @brief AdaDelta optimization process
    InputType optimize() override
    {
        // get params
        auto max_iterations = params->getMaxIterations();
        auto min_gradient = params->getMinGradient();
        auto sgd_times = params->getSGDTimes();
        auto gamma = params->getGamma();
        auto epsilon = params->getEpsilon();

        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/AdaDelta.txt");
        while (true) {
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

                if (params->getIterationTimes() < sgd_times)
                {
                    d = -g;
                    stepsize = ss->search(d);
                    dx = stepsize * d.transpose();

                    for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                    {
                        s[i] = gamma * s_last[i] + (1 - gamma) * g[i] * g[i];
                        s_dx[i] = gamma * s_dx_last[i] + (1 - gamma) * dx[i] *dx[i];
                    }
                }
                else
                {
                    for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                    {
                        s[i] = gamma * s_last[i] + (1 - gamma) * g[i] * g[i];
                        dx[i] = -sqrt(s_dx_last[i]) / (sqrt(s[i]) + epsilon) * g[i];
                        s_dx[i] = gamma * s_dx_last[i] + (1 - gamma) * dx[i] *dx[i];
                    }
                }
                x_next = x + dx;
                // update x
                f->setX(x_next);
                s_last = s;
                s_dx_last = s_dx;

                params->nextIteration();
            }
        }
    }

private:
    AdaDeltaParams* params;
    JacobianType s;         // cumulative sum of squares of gradients
    JacobianType s_last;    // cumulative sum of squares of gradients at last time
    InputType x;            // x at time k
    InputType x_next;       // x at time k+1
    InputType dx;
    InputType s_dx;
    InputType s_dx_last;
    JacobianType g;         // gradient at time k

};
}


#endif
