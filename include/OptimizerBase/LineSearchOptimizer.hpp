#ifndef LINESEARCH_HPP
#define LINESEARCH_HPP

#include <OptimizerBase/OptimizerBase.hpp>

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

    using Base::f;

public:
    virtual ~LineSearchOptimizer()
    {
        //delete stepsize_searcher;
    }

protected:
    //AccurateSearchBase<FunctorType>* stepsize_searcher; // Stepsize search method
    //InaccurateSearchBase<FunctorType>* inaccurate_stepsize;
    StepsizeSearchBase<FunctorType>* ss; // Stepsize searcher
    JacobianType d; //The direction of descent
    T stepsize; // Stepsize for one iteration

};
}

#endif
