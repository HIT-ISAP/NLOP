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
                //this->printProcessInformation();

                x = f->getX();
                g = f->getJacobian();

                v = params->gamma_v * v_last + (1 - params->gamma_v) * g;

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                {
                    s[i] = params->gamma_s * s_last[i] + (1 - params->gamma_s) * g[i] * g[i];
                }

                v_vee = v / (1 - pow(params->gamma_v, params->iteration_times));
                s_vee = s / (1 - pow(params->gamma_s, params->iteration_times));

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                {
                    x_next[i] = x[i] - params->alpha * v_vee[i] / (params->epsilon + sqrt(s_vee[i]));
                }

                f->setX(x_next);

                s_last = s;
                v_last = v;
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
