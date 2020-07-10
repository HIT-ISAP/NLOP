#ifndef CONJUATEGRADIENT_HPP
#define CONJUATEGRADIENT_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/ConjuateGradientParams.hpp>

namespace NLOP {
/// @class ConjuateGradientOptimizer
/// @brief Conjuate gradient method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class ConjuateGradientOptimizer: public LineSearchOptimizer<FunctorType>
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
    /// @brief Constructor and Deconstructor
    ConjuateGradientOptimizer() {}
    ~ConjuateGradientOptimizer() { delete ss; }

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              ConjuateGradientParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        this->initStepsizeMethod(params);
        this->ss->bind(this->f);
        this->ss->setParams(params->stepsize_params);

        beta = 0;
        last_d.setZero(1, InputType::RowsAtCompileTime);
        last_g.setZero(1, InputType::RowsAtCompileTime);
    }

    /// @brief Conjuate Gradient optimization process
    InputType optimize() override
    {
        // get params
        auto min_gradient = params->getMinGradient();
        auto max_iterations = params->getMaxIterations();

        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/conjuate gradient -- "
                              + params->StepsizeMethodTranslator(params->getStepsizeMethod()) + ".txt");
        while (true){
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
                this->printResult(params);
                if (this->writer.is_open())
                    this->writer.close();
                return f->getX();
            }
            else
            {
                g = f->getJacobian();

                // compute beta
                if (params->getIterationTimes() == 0)
                    beta = 0;
                else
                    beta = (g*g.transpose()/(last_g*last_g.transpose()))(0,0);

                // update the direction of descent d
                d = -g + beta*last_d;

                // search stepsize at the direction of d
                stepsize = ss->search(d);

                // update x
                f->setX(f->getX() + stepsize * d.transpose());

                // save gradient and direction at last time
                last_g = g;
                last_d = d;
                params->nextIteration();
            }
        }
        if (this->writer.is_open())
            this->writer.close();
    }

private:
    ConjuateGradientParams* params;
    JacobianType last_g;                // gradient at time k-1
    JacobianType g;                     // gradient at time k
    JacobianType last_d;                // direction at time k-1

    T beta = 0;
};
}

#endif
