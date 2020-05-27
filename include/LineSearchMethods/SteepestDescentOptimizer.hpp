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
              SteepestDescentParams* params, StepsizeSearchBase<FunctorType>* ss)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        this->ss = ss;
        this->ss->init(this->f);
    }

    /// @brief Steepest Descent optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations();
        this->writer.open("../data/"
                          "steepest descent -- WolfePowell.txt");
        while (true) {
            this->updateValueAndJacobian();

            this->writeInformation();
            if (params->iteration_times > params->max_iteration_times)
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
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

                // The direction of descent is negative gradient
                d = -f->getJacobian();

                // Search stepsize at the direction of d
                stepsize = ss->search(d);

                // Update x: x(k+1) = x(k) + stepsize(k) * d(k)
                f->setX(f->getX() + stepsize * d.transpose());
                //this->printProcessInformation();
            }
        }
        this->writer.close();
    }
private:
    SteepestDescentParams* params;
};
}

#endif
