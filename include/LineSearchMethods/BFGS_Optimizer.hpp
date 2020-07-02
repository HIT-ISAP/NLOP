#ifndef BFGS_OPTIMIZER_HPP
#define BFGS_OPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/BFGS_Params.hpp>

namespace NLOP {
/// @class BFGS_Optimizer
/// @brief BFGS Quasi-Newton optimizer
/// @param FunctorType Target funtion type
template<typename FunctorType>
class BFGS_Optimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::HessianType;
    using typename LineSearch::T;

    using LineSearch::f;
    using LineSearch::ss;

    using LineSearch::stepsize;
    using LineSearch::d;

public:
    /// @brief Constructors
    BFGS_Optimizer() {}

    ~BFGS_Optimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              BFGS_Params* params, StepsizeSearchBase<FunctorType>* ss)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        this->ss = ss;
        this->ss->init(this->f);

        delta_x.setZero(InputType::RowsAtCompileTime, 1);
        delta_g.setZero(1, InputType::RowsAtCompileTime);
        H_last.setIdentity(InputType::RowsAtCompileTime, InputType::RowsAtCompileTime);
        I.setIdentity(InputType::RowsAtCompileTime, InputType::RowsAtCompileTime);
    }

    /// @brief BFGS optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations();
        this->writer.open("../data/"
                          "BFGS -- GoldenSection.txt");
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

                // Get gradient
                g = f->getJacobian();

                // Update H
                if (params->iteration_times == 1)
                {
                    H.setIdentity(InputType::RowsAtCompileTime, InputType::RowsAtCompileTime);
                }
                else
                {
                    delta_g = g - g_last;
                    H = (I - delta_x * delta_g / delta_g.dot(delta_x)) * H_last
                            * (I - delta_g.transpose() * delta_x.transpose() / delta_g.dot(delta_x))
                            + delta_x * delta_x.transpose() / delta_g.dot(delta_x);
                }

                // Update direction
                d = -g * H.transpose();

                // Search stepsize along the new direction
                stepsize = ss->search(d);

                // Update x
                delta_x = stepsize * d.transpose();
                f->setX(f->getX() + delta_x);

                // save H, gradient and x at last time
                H_last = H;
                g_last = g;
                x_last = x;

                //this->printProcessInformation();
            }
        }
        this->writer.close();
    }

    BFGS_Params* params; // optimizer params

    JacobianType g; // gradient at time k
    JacobianType g_last; // gradient at time k-1
    JacobianType delta_g; // difference of g between time k and k-1

    InputType x; // x at time k
    InputType x_last; // x at time k-1
    InputType delta_x; // difference of x between time k and k-1

    HessianType H; // H at time k
    HessianType H_last; // H at time k-1

    HessianType I; // Identity matrix
};
}

#endif
