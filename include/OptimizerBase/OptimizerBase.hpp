#ifndef OPTIMIZERBASE_HPP
#define OPTIMIZERBASE_HPP

#include <iostream>
#include <OptimizerParams/OptimizerParamsBase.hpp>
#include <StepsizeSearch/StepsizeSearchBase.hpp>

namespace NLOP {
/// @class NLOP::OptimizerBase
/// @brief Abstract base class for all non-linear optimization methods
/// @param FunctorType Target function type
template<typename FunctorType>
class OptimizerBase
{
protected:
    using T = typename FunctorType::Scalar;
    using InputType = typename FunctorType::InputType;
    using ValueType = typename FunctorType::ValueType;
    using JacobianType = typename FunctorType::JacobianType;
    using HessianType = typename FunctorType::HessianType;

public:
    /// TODO: assert the dimension is correct


    /// @brief Update y with recent x
    void updateValue()
    {
        adjac(f->x, f->y);
    }

    /// @brief Update y and jacobian with recent x
    void updateValueAndJacobian()
    {
        adjac(f->x, f->y, f->jacobian);
    }

    /// @brief Print optimization result
    void printResult()
    {
        std::cout << "Optimization Finished!" << std::endl;
        std::cout << "Optimal x: (" << f->getX().transpose() << ")" << std::endl;
        std::cout << "f(x) = " << f->getY() << std::endl;
        std::cout << "Gradient: (" << f->getJacobian() << ")" << std::endl;
    }

    /// @brief Print optimization process information
    void printProcessInformation()
    {
        std::cout << "x: (" << f->getX().transpose() << ")" << std::endl;
        std::cout << "f(x) = " << f->getY() << std::endl;
        std::cout << "Gradient: (" << f->getJacobian() << ")" << std::endl;
    }

    /// @brief Iteratively compute the optimal x
    virtual InputType optimize() = 0;

    void printInitialConfigurations()
    {
        std::cout << "Initial Configurations: " << "\n"
                  << "x0: (" << f->getX().transpose() << ") \n"
                  << "f(x0) = " << f->getY() << std::endl;
    }

    virtual ~OptimizerBase()
    {
        delete f;
    }

protected:
    FunctorType* f; // Target function
    Eigen::AutoDiffJacobian<FunctorType> adjac; // Tool to compute value and jacobian

};
}

#endif
