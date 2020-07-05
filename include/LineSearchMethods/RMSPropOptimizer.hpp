#ifndef RMSPROPOPTIMIZER_HPP
#define RMSPROPOPTIMIZER_HPP

#include <OptimizerBase/LineSearchOptimizer.hpp>
#include <OptimizerParams/RMSPropParams.hpp>

namespace NLOP {
/// @class RMSPropOptimizer
/// @brief RMSProp method optimizer
/// @param FunctorType Target function type
template<typename FunctorType>
class RMSPropOptimizer: public LineSearchOptimizer<FunctorType>
{
protected:
    using LineSearch = LineSearchOptimizer<FunctorType>;

    using typename LineSearch::InputType;
    using typename LineSearch::ValueType;
    using typename LineSearch::JacobianType;
    using typename LineSearch::T;

    using LineSearch::f;

public:
    RMSPropOptimizer() {}

    /// @brief Initialize
    void init(const InputType& initial, FunctorType* f,
              RMSPropParams* params)
    {
        this->f = f;
        this->f->setX(initial);
        this->updateValue();
        this->params = params;

        s_last.setZero(1, InputType::RowsAtCompileTime);

    }

    /// @brief RMSProp optimization process
    InputType optimize() override
    {
        if (params->getVerbosity() == RMSPropParams::SUMMARY
                 || params->getVerbosity() == RMSPropParams::DETAIL)
        {
            params->print("RMSProp optimization");
            this->printInitialConfigurations();
        }
        //this->writer.open("../data/"
        //                  "RMSProp.txt");
        while (true) {
            this->updateValueAndJacobian();
            //this->writeInformation();
            if (params->getIterationTimes()> params->getMaxIterations())
            {
                std::cerr << "Beyond max iteration times, cannot convergence" << std::endl;
                this->printResult();
                return f->getX();
            }
            if (f->getJacobian().norm() < params->getMinGradient())
            {
                if (params->getVerbosity() == RMSPropParams::SUMMARY
                         || params->getVerbosity() == RMSPropParams::DETAIL)
                {
                    std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
                    this->printResult();
                }
                return f->getX();
            }
            else
            {
                params->nextIteration();

                if (params->getVerbosity() == RMSPropParams::DETAIL)
                {
                    std::cout << "Iteration times: " << params->getIterationTimes() << std::endl;
                    this->printProcessInformation();
                }

                x = f->getX();
                g = f->getJacobian();

                //s = params->gamma * s_last + (1-params->gamma) * (g.dot(g));

                for (int i = 0; i < InputType::RowsAtCompileTime; ++i)
                {
                    s[i] = params->getGamma() * s_last[i] + (1-params->getGamma()) * g[i] * g[i];
                    //s[i] += g[i] * g[i];
                    x_next[i] = x[i] - params->getAlpha() / (params->getEpsilon() + sqrt(s[i])) * g[i];
                }

                f->setX(x_next);
                s_last = s;
            }
        }
        //this->writer.close();
    }

private:
    RMSPropParams* params;
    JacobianType s;         // cumulative sum of squares of gradients
    JacobianType s_last;    // cumulative sum of squares of gradients at last time
    InputType x;            // x at time k
    InputType x_next;       // x at time k+1
    JacobianType g;         // gradient at time k
};
}

#endif
