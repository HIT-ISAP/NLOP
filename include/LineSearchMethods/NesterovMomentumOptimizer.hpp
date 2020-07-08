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
        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/Nesterov momentum.txt");
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
                // compute momentum v
                x_next = f->getX() + params->getBeta() * v_last.transpose();
                computeValueAndJacobian(x_next, &y_next, &jac_next);
                v = params->getBeta() * v_last - params->getAlpha() * jac_next;

                // update x
                f->setX(f->getX() + v.transpose());
                v_last = v;

                params->nextIteration();
            }
        }
        if (this->writer.is_open())
            this->writer.close();
    }

private:
    NesterovMomentumParams* params;
    JacobianType v;         // gradient momentum at time k
    JacobianType v_last;    // gradient momentum at time k-1
    InputType x;

    JacobianType jac_next;  // predicted Jacobian at time k+1
    ValueType y_next;       // predicted f(x) at time k+1
    InputType x_next;       // predicted x at time k+1
};
}

#endif
