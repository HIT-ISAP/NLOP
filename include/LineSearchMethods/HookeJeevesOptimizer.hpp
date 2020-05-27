#ifndef HOOKEJEEVESOPTIMIZER_HPP
#define HOOKEJEEVESOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/HookeJeevesParams.hpp>

namespace NLOP {
/// @class HookeJeevesOptimizer
/// @brief Hooke Jeeves method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class HookeJeevesOptimizer: public LineSearchOptimizer<FunctorType>
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
    HookeJeevesOptimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              HookeJeevesParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        stepsize = params->init_delta;
    }

    /// @brief Hooke Jeeves optimization method
    InputType optimize() override
    {
        this->printInitialConfigurations();
        this->writer.open("../data/"
                          "Hooke Jeeves.txt");
        x = f->getX();
        y = x;
        step1();
        this->writer.close();
    }

private:
    HookeJeevesParams* params;

    InputType x;
    InputType x_next;
    InputType y;
    InputType y_next;

    void step1()
    {
        for (int i = 0; i < InputType::RowsAtCompileTime; i++)
        {
            d.setZero(1, InputType::RowsAtCompileTime);
            d[i] = 1;

            if ((*f)(y + stepsize * d.transpose()) < (*f)(y))
            {
                y_next = y + stepsize * d.transpose();
            }
            else if ((*f)(y - stepsize * d.transpose()) < (*f)(y))
            {
                y_next = y - stepsize * d.transpose();
            }
            else
            {
                y_next = y;
            }

            y = y_next;
        }
        step2();
    }

    void step2()
    {
        if ((*f)(y_next) < (*f)(x))
            step3();
        else
            step4();
    }

    void step3()
    {
        params->iteration_times++;

        x_next = y_next;
        f->setX(x_next);
        this->updateValue();
        this->writeInformation();

        y = x_next + params->alpha * (x_next - x);
        x = x_next;
        step1();
    }

    InputType step4()
    {
        if (stepsize < params->epsilon)
        {
            params->iteration_times++;

            f->setX(x_next);
            this->updateValue();

            std::cout << "Iteration times: " << params->iteration_times << std::endl;
            std::cout << "Optimization Finished!" << std::endl;
            std::cout << "Optimal x: (" << f->getX().transpose() << ")" << std::endl;
            std::cout << "f(x) = " << f->getY() << std::endl;
            return x_next;
        }
        else
        {
            stepsize *= 0.5;
            y = x;
            x_next = x;
            step1();
        }
    }
};
}

#endif
