#ifndef DFP_OPTIMIZER_HPP
#define DFP_OPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/DFP_Params.hpp>

namespace NLOP {
/// @class DFP_Optimizer
/// @brief DFP Quasi-Newton method optimizer
/// @param FunctorType Target funtion type
template<typename FunctorType>
class DFP_Optimizer: public LineSearchOptimizer<FunctorType>
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
    DFP_Optimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              DFP_Params* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        this->initStepsizeMethod(params);
        //this->ss->init(this->f);

        delta_x.setZero(InputType::RowsAtCompileTime, 1);
        delta_g.setZero(1, InputType::RowsAtCompileTime);
        H_last.setIdentity(InputType::RowsAtCompileTime, InputType::RowsAtCompileTime);
    }

    /// @brief DFP optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/DFP -- "
                              + params->StepsizeMethodTranslator(params->getStepsizeMethod()) + ".txt");
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
                // Get gradient
                g = f->getJacobian();

                // Update H
                if (params->getIterationTimes() == 0)
                    H.setIdentity(InputType::RowsAtCompileTime, InputType::RowsAtCompileTime);
                else
                {
                    delta_g = g - g_last;
                    H = H_last + delta_x * delta_x.transpose() / (delta_x.dot(delta_g))
                        - H_last * delta_g.transpose() * delta_g * H_last / ((delta_g * H * delta_g.transpose())(0,0));
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
                params->nextIteration();
            }
        }
        if (this->writer.is_open())
            this->writer.close();
    }

private:
    DFP_Params* params;     // optimizer params

    JacobianType g;         // gradient at time k
    JacobianType g_last;    // gradient at time k-1
    JacobianType delta_g;   // difference of g between time k and k-1

    InputType x;            // x at time k
    InputType x_last;       // x at time k-1
    InputType delta_x;      // difference of x between time k and k-1

    HessianType H;          // H at time k
    HessianType H_last;     // H at time k-1
};
}

#endif
