#ifndef NEWTONOPTIMIZER_HPP
#define NEWTONOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/NewtonParams.hpp>

namespace NLOP {
/// @class NewtonOptimizer
/// @brief Newton method optimizer
/// @param FunctorType Target function type
/// @param HessianFunctorType mammal Hessian functor of target function
template<typename FunctorType, typename HessianFunctorType>
class NewtonOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;

public:
    /// @brief Constructors
    NewtonOptimizer() {}

    /// @brief Initialization
    void init(const InputType& initial, FunctorType* f,
              NewtonParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
    }

    /// @brief Newton optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations();
        while (true){
            this->updateValueAndJacobian();
            if (params->iteration_times > params->max_iteration_times)
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                return f->getX();
            }
            else
            {
                params->iteration_times++;
                //this->printProcessInformation();

                x = f->getX();
                g = f->getJacobian();

                //auto h = H(x);

                delta_x = H(x).inverse() * g.transpose();

                if (delta_x.norm() < params->min_delta_x)
                {
                    std::cout << "Iteration times: " << params->iteration_times << std::endl;
                    this->printResult();
                    return f->getX();
                }

                f->setX(x - delta_x);
            }
        }
    }

private:
    NewtonParams* params;

    HessianFunctorType H;
    InputType delta_x;
    InputType x;
    JacobianType g;
};
}

#endif
