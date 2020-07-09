#ifndef STEEPESTDESCENT_HPP
#define STEEPESTDESCENT_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/SteepestDescentParams.hpp>

namespace NLOP {
/// @class SteepestDescentOptimizer
/// @brief Steepest descent method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class SteepestDescentOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;
    using LineSearch::ss;

    using LineSearch::stepsize;
    using LineSearch::d;

public:
    /// @brief Constructors
    SteepestDescentOptimizer() {}
    ~SteepestDescentOptimizer() {}

    /// @brief Initialize SteepestDescentOptimizer
    void init(const InputType& initial, FunctorType* f,
              SteepestDescentParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        this->initStepsizeMethod(params);
        this->ss->bind(this->f);
        this->ss->setParams(params->stepsize_params);
    }

    /// @brief Steepest Descent optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/steepest descent -- "
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
                // The direction of descent is negative gradient
                d = -f->getJacobian();

                // Search stepsize at the direction of d
                stepsize = ss->search(d);

                // Update x: x(k+1) = x(k) + stepsize(k) * d(k)
                f->setX(f->getX() + stepsize * d.transpose());
                params->nextIteration();
            }
        }
        if (this->writer.is_open())
            this->writer.close();
    }
private:
    SteepestDescentParams* params;
};
}

#endif
