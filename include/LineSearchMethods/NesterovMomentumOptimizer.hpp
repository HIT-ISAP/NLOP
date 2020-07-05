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
        if (params->getVerbosity() == NesterovMomentumParams::SUMMARY
                 || params->getVerbosity() == NesterovMomentumParams::DETAIL)
        {
            params->print("Nesterov Momentum optimization");
            this->printInitialConfigurations();
        }
        //this->writer.open("../data/"
        //                  "nesterov momentum.txt");
        while (true) {
            this->updateValueAndJacobian();
            //this->writeInformation();
            if (params->getIterationTimes() > params->getMaxIterations())
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                this->printResult();
                return f->getX();
            }
            if (f->getJacobian().norm() < params->getMinGradient())
            {
                if (params->getVerbosity() == NesterovMomentumParams::SUMMARY
                         || params->getVerbosity() == NesterovMomentumParams::DETAIL)
                {
                    std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
                    this->printResult();
                }
                return f->getX();
            }
            else
            {
                params->nextIteration();

                // compute momentum v
                x_next = f->getX() + params->getBeta() * v_last.transpose();
                computeValueAndJacobian(x_next, &y_next, &jac_next);
                v = params->getBeta() * v_last - params->getAlpha() * jac_next;

                f->setX(f->getX() + v.transpose());
                v_last = v;

                if (params->getVerbosity() == NesterovMomentumParams::DETAIL)
                {
                    std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
                    this->printProcessInformation();
                }
            }
        }
        //this->writer.close();
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
