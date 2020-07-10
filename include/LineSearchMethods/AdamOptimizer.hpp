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
        // get params
        auto max_iterations = params->getMaxIterations();
        auto min_gradient = params->getMinGradient();
        auto alpha = params->getAlpha();
        auto gamma_v = params->getGammaV();
        auto gamma_s = params->getGammaS();
        auto epsilon = params->getEpsilon();

        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/Adam.txt");
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

                v = gamma_v * v_last + (1 - gamma_v) * g;

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                    s[i] = gamma_s * s_last[i] + (1 - gamma_s) * g[i] * g[i];

                v_vee = v / (1 - pow(gamma_v, params->getIterationTimes()+1));
                s_vee = s / (1 - pow(gamma_s, params->getIterationTimes()+1));

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                    x_next[i] = x[i] - alpha * v_vee[i] / (epsilon + sqrt(s_vee[i]));

                // update x
                f->setX(x_next);

                s_last = s;
                v_last = v;
                params->nextIteration();
            }
        }
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
