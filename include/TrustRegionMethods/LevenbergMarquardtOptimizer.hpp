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
        this->printInitialConfigurations(params);
        if (params->isLogFile())
            this->writer.open("../data/LevenbergMarquardt.txt");
        while (true) {
            this->updateValueAndJacobian();
            this->printProcessInformation(params);
            if (this->writer.is_open())
                this->writeInformation();
            if (params->getIterationTimes() > params->getMaxIterations())
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                return f->getX();
            }
            else
            {
                x = f->getX();
                g = f->getJacobian();

                H_m = H(x) + epsilon * I;

                if (!isPosDef()) // if H_m is not positive definite, increase damping factor
                {
                    epsilon *= 4;
                    continue;
                }

                delta_x = H_m.llt().solve(-g.transpose());

                if (delta_x.norm() < params->getMinDeltaX())
                {
                    this->printResult(params);
                    return f->getX();
                }

                x_next = x + delta_x;

                // compute the error of Taylor expansion
                R = ((*f)(x_next) - (*f)(x)) /
                        (0.5 * delta_x.transpose() * H(x) * delta_x + g * delta_x)(0, 0);

                if (R < 0)
                    ; /// @todo
                else if (R < 0.25)   // the error of Taylor expansion is too big, reduce the stepsize
                    epsilon *= 4;
                else if (R > 0.75)   // the error of Taylor expansion is acceptable, increase the stepsize
                    epsilon *= 0.5;

                // update x
                f->setX(x_next);
                params->nextIteration();
            }
        }
        if (this->writer.is_open())
            this->writer.close();
    }

private:
    LevenbergMarquardtParams* params;

    T epsilon;              // damping factor
    T R;                    // param to describe the fitting degree of second order Taylor expansion
    HessianFunctorType H;   // Hessian matrix
    HessianType H_m;        // damped hessian matrix
    HessianType I;          // identity matrix
    InputType delta_x;      // delta x from time k to k+1
    InputType x;            // x at time k
    InputType x_next;       // x at time k+1
    JacobianType g;         // Jacobian

    // check a square matrix whether is positive definite
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
