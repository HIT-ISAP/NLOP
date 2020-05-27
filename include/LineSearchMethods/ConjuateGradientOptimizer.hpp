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
    /// @brief Constructors
    ConjuateGradientOptimizer() {}

    ~ConjuateGradientOptimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              ConjuateGradientParams* params, StepsizeSearchBase<FunctorType>* ss)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        this->ss = ss;
        this->ss->init(this->f);

        beta = 0;
        last_d.setZero(1, InputType::RowsAtCompileTime);
        last_g.setZero(1, InputType::RowsAtCompileTime);
    }

    /// @brief Conjuate Gradient optimization process
    InputType optimize() override
    {
        this->printInitialConfigurations();
        this->writer.open("../data/"
                          "conjuate gradient -- WolfePowell.txt");
        while (true){
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

                // Compute beta
                if (params->iteration_times == 1)
                    beta = 0;
                else
                    beta = (g*g.transpose()/(last_g*last_g.transpose()))(0,0);

                // Update the direction of descent d
                d = -g + beta*last_d;

                // Search stepsize at the direction of d
                stepsize = ss->search(d);

                f->setX(f->getX() + stepsize * d.transpose());

                // save gradient and direction at last time
                last_g = g;
                last_d = d;

                //this->printProcessInformation();
            }
        }
        this->writer.close();
    }

    ConjuateGradientParams* params;
    JacobianType last_g ; // gradient at time k-1
    JacobianType g; // gradient at time k
    JacobianType last_d; // direction at time k-1

    T beta = 0;
};
}

#endif
