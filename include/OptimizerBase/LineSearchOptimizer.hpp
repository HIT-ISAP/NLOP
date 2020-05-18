#ifndef LINESEARCH_HPP
#define LINESEARCH_HPP

#include <OptimizerBase/OptimizerBase.hpp>
#include <StepsizeSearch/Accurate/DichotomousMethod.hpp>
#include <StepsizeSearch/Accurate/FibonacciMethod.hpp>
#include <StepsizeSearch/Accurate/GoldenSectionMethod.hpp>
#include <StepsizeSearch/Inaccurate/GoldsteinMethod.hpp>

namespace NLOP {

/// @class NLOP::LineSearchOptimizer
/// @brief Abstract base class for all line search methods
/// @param T The numeric scalar type
/// @param N The dimension of variable x
template<typename T, int N, typename FunctorType, typename PhiFunctortype>
class LineSearchOptimizer: public OptimizerBase<T, N, FunctorType>
{
protected:
    using Base = OptimizerBase<T, N, FunctorType>;
    using Base::f;
    using Base::params;
    using typename Base::InputType;

public:
    virtual ~LineSearchOptimizer()
    {
        delete stepsize_searcher;
    }

protected:
    AccurateSearchBase<T, PhiFunctortype>* stepsize_searcher; // Stepsize search method
    InaccurateSearchBase<FunctorType>* inaccurate_stepsize;

    T stepsize; // Stepsize for one iteration

};
}

#endif
