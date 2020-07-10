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
    /// @brief Constructor and Deconstructor
    BFGS_Optimizer() {}
    ~BFGS_Optimizer() { delete ss; }

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              BFGS_Params* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        this->initStepsizeMethod(params);
        this->ss->bind(this->f);
        this->ss->setParams(params->stepsize_params);

        delta_x.setZero(InputType::RowsAtCompileTime, 1);
        delta_g.setZero(1, InputType::RowsAtCompileTime);
        H_last.setIdentity(InputType::RowsAtCompileTime, InputType::RowsAtCompileTime);
        I.setIdentity(InputType::RowsAtCompileTime, InputType::RowsAtCompileTime);
    }

    /// @brief BFGS optimization process
    InputType optimize() override
    {
        // get params
        auto max_iterations = params->getMaxIterations();
        auto min_gradient = params->getMinGradient();

        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/BFGS -- "
                              + params->StepsizeMethodTranslator(params->getStepsizeMethod()) + ".txt");
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
                if (this->writer.is_open())
                    this->writer.close();
                this->printResult(params);
                return f->getX();
            }
            else
            {
                // get gradient
                g = f->getJacobian();

                // update H
                if (params->getIterationTimes() == 0)
                    H.setIdentity(InputType::RowsAtCompileTime, InputType::RowsAtCompileTime);
                else
                {
                    delta_g = g - g_last;
                    H = (I - delta_x * delta_g / delta_g.dot(delta_x)) * H_last
                            * (I - delta_g.transpose() * delta_x.transpose() / delta_g.dot(delta_x))
                            + delta_x * delta_x.transpose() / delta_g.dot(delta_x);
                }
                // update direction
                d = -g * H.transpose();

                // search stepsize along the new direction
                stepsize = ss->search(d);

                // update x
                delta_x = stepsize * d.transpose();
                f->setX(f->getX() + delta_x);

                // save H, gradient and x at last time
                H_last = H;
                g_last = g;
                x_last = x;

                params->nextIteration();
            }
        }
    }

private:
    BFGS_Params* params;    // optimizer params

    JacobianType g;         // gradient at time k
    JacobianType g_last;    // gradient at time k-1
    JacobianType delta_g;   // difference of g between time k and k-1

    InputType x;            // x at time k
    InputType x_last;       // x at time k-1
    InputType delta_x;      // difference of x between time k and k-1

    HessianType H;          // H at time k
    HessianType H_last;     // H at time k-1

    HessianType I;          // identity matrix
};
}

#endif
