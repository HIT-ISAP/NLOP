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
        if (params->getVerbosity() == NewtonParams::SUMMARY
                 || params->getVerbosity() == NewtonParams::DETAIL)
        {
            params->print("Newton's method optimization");
            this->printInitialConfigurations();
        }
        this->writer.open("../data/"
                          "Newton.txt");
        while (true){
            this->updateValueAndJacobian();
            this->writeInformation();
            if (params->getIterationTimes() > params->getMaxIterations())
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                return f->getX();
            }
            else
            {
                params->nextIteration();
                if (params->getVerbosity() == NewtonParams::DETAIL)
                {
                    std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
                    this->printProcessInformation();
                }

                x = f->getX();
                g = f->getJacobian();

                delta_x = H(x).inverse() * g.transpose();

                if (delta_x.norm() < params->getMinDeltaX())
                {
                    if (params->getVerbosity() == NewtonParams::SUMMARY
                             || params->getVerbosity() == NewtonParams::DETAIL)
                    {
                        std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
                        this->printResult();
                    }
                    return f->getX();
                }
                f->setX(x - delta_x);
            }
        }
        this->writer.close();
    }

private:
    NewtonParams* params;

    HessianFunctorType H; // Hessian matrix
    InputType delta_x;    // Delta x for every step
    InputType x;
    JacobianType g;       // Jacobian
};
}

#endif
