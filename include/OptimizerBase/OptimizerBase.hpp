#ifndef OPTIMIZERBASE_HPP
#define OPTIMIZERBASE_HPP

#include <iostream>
#include <fstream>
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
    void printResult(OptimizerParamsBase* params)
    {
        if (params->getVerbosity() == OptimizerParamsBase::SUMMARY
                 || params->getVerbosity() == OptimizerParamsBase::DETAIL)
        {
            std::cout << "Optimization Finished!" << "\n";
            std::cout << "Iteration times: " << params->getIterationTimes() << "\n";
            std::cout << "Optimal x: (" << f->getX().transpose() << ")" << "\n";
            std::cout << "f(x) = " << f->getY() << "\n";
            std::cout << "Gradient: (" << f->getJacobian() << ")" << std::endl;
        }
    }

    /// @brief Print optimization process information
    void printProcessInformation(OptimizerParamsBase* params)
    {
        if (params->getVerbosity() == OptimizerParamsBase::DETAIL)
        {
            std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
            std::cout << "x: (" << f->getX().transpose() << ")" << "\n";
            std::cout << "f(x) = " << f->getY() << "\n";
            std::cout << "Gradient: (" << f->getJacobian() << ")" << "\n";
            std::cout << "*********************************************" << std::endl;
        }
    }

    /// @brief Print optimizer initial configuration
    void printInitialConfigurations(OptimizerParamsBase* params)
    {
        if (params->getVerbosity() == OptimizerParamsBase::SUMMARY
                 || params->getVerbosity() == OptimizerParamsBase::DETAIL)
        {
            params->print();
            std::cout << "Initial Configurations: " << "\n"
                      << "x0: (" << f->getX().transpose() << ") \n"
                      << "f(x0) = " << f->getY() << "\n"
                      << "*********************************************" << std::endl;
        }
    }

    /// @brief Write information into log files
    void writeInformation()
    {
        writer << f->getX()[0] << " "
               << f->getX()[1] << " "
               << f->getY() << "\n";
    }

    /// @brief Iteratively compute the optimal x
    virtual InputType optimize() = 0;

    virtual ~OptimizerBase() { delete f; }

protected:
    FunctorType* f;                             // Target function
    Eigen::AutoDiffJacobian<FunctorType> adjac; // Tool to compute value and jacobian
    std::ofstream writer;                       // Tool to write information into txt

};
}

#endif
