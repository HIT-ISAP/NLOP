#ifndef LM_OPTIMIZER_HPP
#define LM_OPTIMIZER_HPP

#include <OptimizerBase/TrustRegionOptimizer.hpp>
#include <OptimizerParams/LevenbergMarquardtParams.hpp>

namespace NLOP {
/// @class LevenbergMarquardtOptimizer
/// @brief Levenberg-Marquardt method optimizer
/// @param FunctorType Target function type
/// @param HessianFunctorType mammal Hessian functor of target function
template<typename FunctorType, typename HessianFunctorType>
class LevenbergMarquardtOptimizer: public TrustRegionOptimizer<FunctorType>
{
protected:
    using TrustRegion = TrustRegionOptimizer<FunctorType>;

    using typename TrustRegion::InputType;
    using typename TrustRegion::ValueType;
    using typename TrustRegion::JacobianType;
    using typename TrustRegion::HessianType;
    using typename TrustRegion::T;

    using TrustRegion::f;

public:
    /// @brief Constructor
    LevenbergMarquardtOptimizer() {}

    /// @brief Initialization
    void init(const InputType& initial, FunctorType* f,
              LevenbergMarquardtParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;
        epsilon = params->getInitEpsilon();
        I.setIdentity(HessianType::RowsAtCompileTime, HessianType::ColsAtCompileTime);
    }

    /// @brief L-M optimization process
    InputType optimize() override
    {
        if (params->getVerbosity() == LevenbergMarquardtParams::SUMMARY
                 || params->getVerbosity() == LevenbergMarquardtParams::DETAIL)
        {
            params->print("L-M optimization");
            this->printInitialConfigurations();
        }
        this->writer.open("../data/"
                          "LevenbergMarquardt.txt");
        while (true) {
            this->updateValueAndJacobian();
            this->writeInformation();
            if (params->getIterationTimes() > params->getMaxIterations())
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                return f->getX();
            }
            else
            {
                params->nextIteration();
                if (params->getVerbosity() == LevenbergMarquardtParams::DETAIL)
                {
                    std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
                    this->printProcessInformation();
                }

                x = f->getX();
                g = f->getJacobian();

                H_m = H(x) + epsilon * I;

                if (!isPosDef())
                {
                    epsilon *= 4;
                    continue;
                }

                delta_x = H_m.llt().solve(-g.transpose());

                if (delta_x.norm() < params->getMinDeltaX())
                {
                    if (params->getVerbosity() == LevenbergMarquardtParams::SUMMARY
                             || params->getVerbosity() == LevenbergMarquardtParams::DETAIL)
                    {
                        std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
                        this->printResult();
                    }
                    return f->getX();
                }

                x_next = x + delta_x;

                R = ((*f)(x_next) - (*f)(x)) /
                        (0.5 * delta_x.transpose() * H(x) * delta_x + g * delta_x)(0, 0);

                if (R < 0)
                    ; /// @todo
                else if (R < 0.25)
                    epsilon *= 4;
                else if (R > 0.75)
                    epsilon *= 0.5;

                f->setX(x_next);
            }
        }
        this->writer.close();
    }

private:
    LevenbergMarquardtParams* params;

    T epsilon;
    T R;

    HessianFunctorType H;

    HessianType H_m;
    HessianType I;

    InputType delta_x;
    InputType x;
    InputType x_next;

    JacobianType g;

    bool isPosDef()
    {
        Eigen::EigenSolver<Eigen::MatrixXd> es(H_m);
        Eigen::VectorXd lambda = es.eigenvalues().real();
        auto min_lambda = lambda.minCoeff();
        if (min_lambda > 0)
            return true;
        else
            return false;
    }
};

}

#endif
