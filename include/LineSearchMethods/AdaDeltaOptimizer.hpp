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

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              AdaDeltaParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        ss = new GoldSectionMethod<FunctorType>;
        ss->init(this->f);

        s.setZero(1, InputType::RowsAtCompileTime);
        s_last.setZero(1, InputType::RowsAtCompileTime);
        s_dx.setZero(InputType::RowsAtCompileTime, 1);
        s_dx_last.setZero(InputType::RowsAtCompileTime, 1);
    }

    /// @brief AdaDelta optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/AdaDelta.txt");
        while (true) {
            this->updateValueAndJacobian();
            this->printProcessInformation(params);
            if (this->writer.is_open())
                this->writeInformation();
            if (params->getIterationTimes() > params->getMaxIterations())
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                return f->getX();
            }
            if (f->getJacobian().norm() < params->getMinGradient())
            {
                this->printResult(params);
                return f->getX();
            }
            else
            {
                x = f->getX();
                g = f->getJacobian();

                if (params->getIterationTimes() < params->getSGDTimes())
                {
                    d = -g;
                    stepsize = ss->search(d);
                    dx = stepsize * d.transpose();

                    for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                    {
                        s[i] = params->getGamma() * s_last[i] + (1 - params->getGamma()) * g[i] * g[i];
                        s_dx[i] = params->getGamma() * s_dx_last[i] + (1 - params->getGamma()) * dx[i] *dx[i];
                    }
                }
                else
                {
                    for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                    {
                        s[i] = params->getGamma() * s_last[i] + (1 - params->getGamma()) * g[i] * g[i];
                        dx[i] = -sqrt(s_dx_last[i]) / (sqrt(s[i]) + params->getEpsilon()) * g[i];
                        s_dx[i] = params->getGamma() * s_dx_last[i] + (1 - params->getGamma()) * dx[i] *dx[i];
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
        if (this->writer.is_open())
            this->writer.close();
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
