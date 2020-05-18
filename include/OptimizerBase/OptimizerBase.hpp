#ifndef OPTIMIZERBASE_HPP
#define OPTIMIZERBASE_HPP

#include <iostream>

#include <OptimizerParams/OptimizerParamsBase.hpp>

namespace NLOP {

/// @class NLOP::OptimizerBase
/// @brief Abstract base class for all non-linear optimization methods
/// @param T The numeric scalar type
/// @param N The dimension of variable x
template<typename T, int N, typename FunctorType>
class OptimizerBase
{
public:
    using InputType = typename Eigen::Matrix<T, N, 1>;
    using ValueType = typename FunctorType::ValueType;
public:
    /// TODO: assert the dimension is correct

    /// Constructor
    OptimizerBase()
    {
        f = new FunctorType;
        params = new OptimizerParamsBase;
    }

    /// @brief Initialize target function f(), variable x
    ///        Choose one of the stepsize search method (default = GoldenSection)
    void init(const InputType& initial, FunctorType* f,
                      OptimizerParamsBase* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->params = params;
    }

    void updateValue()
    {
        adjac(f->x, f->y);
    }

    void updateValueAndJacobian()
    {
        adjac(f->x, f->y, f->jacobian);
    }

    /// @brief Print optimization result
    void printResult()
    {
        std::cout << "Gradient: (" << f->getJacobian() << ")" << std::endl;
        std::cout << "Steepest Descent Optimization Finished!" << std::endl;
        std::cout << "Optimal x: (" << f->getX().transpose() << ")" << std::endl;
        std::cout << "f(x) = " << f->getY() << std::endl;
    }

    /// @brief Print optimization process information
    void printProcessInformation()
    {
        std::cout << "Gradient: " << "\n" << f->getJacobian() << std::endl;
        std::cout << "x: " << "\n" << f->getX() << std::endl;
        std::cout << "f(x) = " << f->getY() << std::endl;
    }

    /// @brief Iteratively compute the optimal x
    virtual InputType optimize() = 0;

    virtual ~OptimizerBase()
    {
        delete f;
        delete params;
    }

protected:
    FunctorType* f; // Target function

    OptimizerParamsBase* params; // Optimizer parameters

    // Tool to compute value and jacobian
    Eigen::AutoDiffJacobian<FunctorType> adjac;

};
}

#endif
