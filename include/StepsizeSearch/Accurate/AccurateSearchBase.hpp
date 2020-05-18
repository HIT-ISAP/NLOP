#ifndef ACCURATESEARCHBASE_HPP
#define ACCURATESEARCHBASE_HPP

#include <StepsizeSearch/StepsizeSearchBase.hpp>

namespace NLOP{

/// @class NLOP::AccurateSearchBase
/// @brief Abstract base class for all accurate search methods
/// @param FunctorType Target function
template<typename FunctorType>
class AccurateSearchBase: public StepsizeSearchBase<FunctorType>
{
protected:
    using Base = StepsizeSearchBase<FunctorType>;
    using typename Base::T;
    using typename Base::InputType;
    using typename Base::ValueType;
    using typename Base::JacobianType;

    using Base::max_iteration_times;
    using Base::iteration_times;
    using Base::f;

public:
    /// @brief set params for one dimensional search methods
    /// @param init_interval 2-dims Vector, the initial interval to search in
    /// @param epsilon The threshold to stop the iteration
    /// @param max_iteration_times The threshold of most iteration times
    ///        (Except Dichotomous Method need to manual set, default = 1000)
    void init(FunctorType* f, Vector<T, 2> init_interval, T epsilon_value)
    {
        this->f = f;
        alpha = init_interval[0];
        beta = init_interval[1];
        epsilon = epsilon_value;
        iteration_times = 0;
    }

    /// @brief Print process information
    virtual void printProcess() = 0;


    /// @brief Initial with default params
    void init(FunctorType* f) override
    {
        this->f = f;
        alpha = 0;
        beta = 1;
        epsilon = 0.001;
        iteration_times = 0;
    }

    /// @brief Reset observation variables
    virtual void reset()
    {
        alpha = 0;
        beta = 1;
        iteration_times = 0;
    }

    /// @brief print one dimensional search information
    //virtual void print();

    /// @brief get current interval [alpha, beta]
    Vector<T, 2> getCurrentInterval()
    {
        Vector<T, 2> current_interval(alpha, beta);
        return current_interval;
    }

    // Observation variables
    T alpha;
    T beta;
    T lambda;

    T epsilon; // Threshold to stop the iteration

};

}

#endif
