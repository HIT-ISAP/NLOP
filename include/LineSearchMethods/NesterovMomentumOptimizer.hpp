#ifndef NESTEROVOPTIMIZER_HPP
#define NESTEROVOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/NesterovMomentumParams.hpp>

namespace NLOP {
/// @class NLOP::NesterovMomentumOptimizer
/// @brief Params of Nestero momentum optimization methods
template<typename FunctorType>
class NesterovMomentumOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;
    using LineSearch::Base::adjac;

public:
    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              NesterovMomentumParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;

        v_last.setZero(1, InputType::RowsAtCompileTime);
    }

    /// @brief Compute y and jacobian with given x
    void computeValueAndJacobian(const InputType& x, ValueType* v, JacobianType* jac)
    {
        adjac(x, v, jac);
    }

    /// @brief Nesterov momentum optimization process
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

                // compute momentum v
                x_next = f->getX() + params->beta * v_last.transpose();
                computeValueAndJacobian(x_next, &y_next, &jac_next);
                v = params->beta * v_last - params->alpha * jac_next;

                f->setX(f->getX() + v.transpose());
                v_last = v;

                //this->printProcessInformation();
            }
        }
    }

private:
    NesterovMomentumParams* params;
    JacobianType v; // gradient momentum at time k
    JacobianType v_last; // gradient momentum at time k-1
    InputType x;
    //JacobianType v_next;

    JacobianType jac_next;
    ValueType y_next;
    InputType x_next;
};
}

#endif
