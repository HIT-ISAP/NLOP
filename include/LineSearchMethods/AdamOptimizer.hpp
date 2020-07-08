#ifndef ADAMOPTIMIZER_HPP
#define ADAMOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/AdamParams.hpp>

namespace NLOP {
/// @class NLOP::AdamOptimizer
/// @brief Adam method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class AdamOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;

public:
    AdamOptimizer() {}
    ~AdamOptimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              AdamParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;

        s_last.setZero(1, InputType::RowsAtCompileTime);
        v_last.setZero(1, InputType::RowsAtCompileTime);
    }

    /// @brief Adam optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/Adam.txt");
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

                v = params->getGammaV() * v_last + (1 - params->getGammaV()) * g;

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                {
                    s[i] = params->getGammaS() * s_last[i] + (1 - params->getGammaS()) * g[i] * g[i];
                }

                v_vee = v / (1 - pow(params->getGammaV(), params->getIterationTimes()+1));
                s_vee = s / (1 - pow(params->getGammaS(), params->getIterationTimes()+1));

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                {
                    x_next[i] = x[i] - params->getAlpha() * v_vee[i] / (params->getEpsilon() + sqrt(s_vee[i]));
                }

                // update x
                f->setX(x_next);

                s_last = s;
                v_last = v;

                params->nextIteration();
            }
        }
        if (this->writer.is_open())
            this->writer.close();
    }

private:
    AdamParams* params;
    JacobianType s;
    JacobianType s_last;
    JacobianType s_vee;

    JacobianType v;
    JacobianType v_last;
    JacobianType v_vee;

    InputType x;
    InputType x_next;

    JacobianType g;

};
}


#endif
