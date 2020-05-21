#ifndef NEWTONOPTIMIZER_HPP
#define NEWTONOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/NewtonParams.hpp>

namespace NLOP {
/// @class NewtonOptimizer
/// @brief Newton method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class NewtonOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;
    using typename LineSearch::Base::HessianType;

    using LineSearch::f;

public:
    /// @brief Constructors
    NewtonOptimizer() {}

    void init(const InputType& initial, FunctorType* f,
              NewtonParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
    }

private:
    NewtonParams* params;

    HessianType H;

    void updateHessian()
    {
        adjac(f->x, f->)
    }
};
}

#endif
