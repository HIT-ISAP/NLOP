#ifndef LINESEARCH_HPP
#define LINESEARCH_HPP

#include <OptimizerBase/OptimizerBase.hpp>

#include <OptimizerParams/LineSearchParams.hpp>

#include <StepsizeSearch/Accurate/DichotomousMethod.hpp>
#include <StepsizeSearch/Accurate/FibonacciMethod.hpp>
#include <StepsizeSearch/Accurate/GoldenSectionMethod.hpp>

#include <StepsizeSearch/Inaccurate/GoldsteinMethod.hpp>
#include <StepsizeSearch/Inaccurate/ArmijoMethod.hpp>
#include <StepsizeSearch/Inaccurate/WolfePowellMethod.hpp>

namespace NLOP {
/// @class NLOP::LineSearchOptimizer
/// @brief Abstract base class for all line search methods
/// @param FunctorType Target function type
template<typename FunctorType>
class LineSearchOptimizer: public OptimizerBase<FunctorType>
{
protected:
    using Base = OptimizerBase<FunctorType>;

    using typename Base::T;

    using typename Base::InputType;
    using typename Base::ValueType;
    using typename Base::JacobianType;
    using typename Base::HessianType;

    using Base::f;

public:
    ~LineSearchOptimizer() {}

    /// @brief Initialize stepsize search method
    void initStepsizeMethod(LineSearchParams* params)
    {
        switch (params->getStepsizeMethod()) {
        case LineSearchParams::GOLDENSECTION:
            ss = new GoldSectionMethod<FunctorType>;
            break;
        case LineSearchParams::DICHOTOMOUS:
            ss = new DichotomousMethod<FunctorType>;
            break;
        case LineSearchParams::FIBONACCI:
            ss = new FibonacciMethod<FunctorType>;
            break;
        case LineSearchParams::ARMIJO:
            ss = new ArmijoMethod<FunctorType>;
            break;
        case LineSearchParams::GOLDSTEIN:
            ss = new GoldsteinMethod<FunctorType>;
            break;
        case LineSearchParams::WOLFEPOWELL:
            ss = new WolfePowellMethod<FunctorType>;
            break;
        default:
            ss = new WolfePowellMethod<FunctorType>;
            break;
        }
    }

protected:
    StepsizeSearchBase<FunctorType>* ss;    // stepsize searcher
    JacobianType d;                         // the direction of descent
    T stepsize;                             // stepsize for one iteration

};
}

#endif
