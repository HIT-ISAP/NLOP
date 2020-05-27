#ifndef MOMENTUMOPTIMIZER_HPP
#define MOMENTUMOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/MomentumParams.hpp>

namespace NLOP {
/// @class MomentumOptimizer
/// @brief Momentum method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class MomentumOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;

public:
    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              MomentumParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;

        v_last.setZero(1, InputType::RowsAtCompileTime);
    }

    /// @brief Momentum optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations();
        this->writer.open("../data/"
                          "Momentum.txt");
        while (true) {
            this->updateValueAndJacobian();
            this->writeInformation();
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

                g = f->getJacobian();

                // compute momentum v
                v = params->beta * v_last - params->alpha * g;

                f->setX(f->getX() + v.transpose());
                v_last = v;

                //this->printProcessInformation();
            }
        }
        this->writer.close();
    }

private:
    MomentumParams* params;
    JacobianType v; // gradient momentum at time k
    JacobianType g; // gradient at time k
    JacobianType v_last; // gradient momentum at time k-1

};
}

#endif
