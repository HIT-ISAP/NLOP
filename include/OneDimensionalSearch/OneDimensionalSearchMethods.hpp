#ifndef ONEDIMENSIONALSEARCHMETHODS_HPP
#define ONEDIMENSIONALSEARCHMETHODS_HPP

#include "Utils/Functor.hpp"
#include "Utils/Matrix.hpp"
#include <iostream>

namespace NLOP{

/// @class NLOP::OneDimSearch
/// @brief Abstract base class for all one dimensional search methods
/// @param T The numeric scalar type
template<typename T, typename FunctorType>
class OneDimSearch
{
public:

    /// @brief Iteratively search the optimal stepsize
    virtual T search() = 0;

    /// @brief set params for one dimensional search methods
    ///
    /// @param init_interval 2-dims Vector, the initial interval to search in
    /// @param epsilon The threshold to stop the iteration
    /// @param max_iteration_times The threshold of most iteration times
    ///        (Except Dichotomous Method need to manual set, default = 1000)
    void init(Vector<T, 2> init_interval, T epsilon_value)
    {
        alpha = init_interval[0];
        beta = init_interval[1];
        epsilon = epsilon_value;
        iteration_times = 0;
    }

    void init()
    {
        alpha = 0;
        beta = 1;
        epsilon = 0.001;
        iteration_times = 0;
    }

    void setPhi(const FunctorType& phi_functor)
    {
        phi = phi_functor;
    }

    /// @brief Set max iteration times
    void setMaxIterations(int value)
    {
        max_iteration_times = value;
    }

    /// @brief print one dimensional search information
    //virtual void print();

    /// @brief get current interval [alpha, beta]
    Vector<T, 2> getCurrentInterval()
    {
        Vector<T, 2> current_interval(alpha, beta);
        return current_interval;
    }

    /// @brief Constructors
    OneDimSearch() {}
    //virtual ~OneDimSearch(){}
    //OneDimSearch(Vector<T, 2> init_interval): alpha(init_interval[0]), beta(init_interval[1]) {}


    /// phi(lambda)
    FunctorType phi;

    T alpha;
    T beta;
    T lambda;
    T epsilon;

    int iteration_times = 0;
    int max_iteration_times = 15;

};

}

#endif
