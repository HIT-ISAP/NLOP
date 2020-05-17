#ifndef LINESEARCH_HPP
#define LINESEARCH_HPP

#include <OptimizerBase/OptimizerBase.hpp>

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
    OneDimSearch<T, PhiFunctortype>* stepsize_searcher; // Stepsize search method

    T stepsize; // Stepsize for one iteration

};
}

#endif
