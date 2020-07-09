#include <iostream>
#include <Utils/Functor.hpp>
#include <LineSearchMethods/SteepestDescentOptimizer.hpp>
#include <LineSearchMethods/ConjuateGradientOptimizer.hpp>
#include <LineSearchMethods/MomentumOptimizer.hpp>
#include <LineSearchMethods/NesterovMomentumOptimizer.hpp>
#include <LineSearchMethods/AdagradOptimizer.hpp>
#include <LineSearchMethods/RMSPropOptimizer.hpp>
#include <LineSearchMethods/AdaDeltaOptimizer.hpp>
#include <LineSearchMethods/AdamOptimizer.hpp>
#include <LineSearchMethods/HookeJeevesOptimizer.hpp>
#include <LineSearchMethods/NewtonOptimizer.hpp>
#include <TrustRegionMethods/LevenbergMarquardtOptimizer.hpp>
#include <LineSearchMethods/DFP_Optimizer.hpp>
#include <LineSearchMethods/BFGS_Optimizer.hpp>

using namespace std;

struct RosenbrockFunctor: public NLOP::Functor<double, 2>
{
    void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const
    {
        (*v)[0] = (1 - x[0])*(1 - x[0]) + 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]);
    }

    void operator ()(const ActiveInput& x, ActiveValue* v) const
    {
        (*v)[0] = (1 - x[0])*(1 - x[0]) + 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]);
    }

    double operator ()(const InputType& x)
    {
        double result = (1 - x[0])*(1 - x[0]) + 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]);
        return result;
    }
};

struct RosenbrockHessianFunctor: public NLOP::Functor<double, 2>
{
    RosenbrockHessianFunctor() {}

    HessianType operator() (const InputType& x)
    {
        HessianType H;
        H(0, 0) = 1200 * x[0] * x[0] - 400 * x[1] + 2;
        H(0, 1) = -400 * x[0];
        H(1, 0) = -400 * x[0];
        H(1, 1) = 200;
        return H;
    }
};

using namespace NLOP;

int main()
{
    RosenbrockFunctor* f = new RosenbrockFunctor;
    Eigen::Vector2d initial_x(0, 0);

    //SteepestDescentOptimizer<RosenbrockFunctor> optimizer;
    //SteepestDescentParams* params= new SteepestDescentParams;
    //params->setVerbosity("DETAIL");
    //params->setMaxItertaions(1000);
    //params->setMinGradient(0.1);
    //params->setStepsizeMethod("GOLDSTEIN");
    //params->setStepsizeAccuracy(0.001);
    //params->setStepsizeIncreaseFactor(3);
    //params->setStepsizeDecreaseFactor(0.3);
    //params->setStepsizeUpperBound(0.01);
    //params->setStepsizeLowerBound(0);
    //params->enableLog();
    //optimizer.init(initial_x, f, params);

    //ConjuateGradientOptimizer<RosenbrockFunctor> optimizer;
    //ConjuateGradientParams* params = new ConjuateGradientParams;
    //params->setStepsizeMethod("GOLDSTEIN");
    //params->setStepsizeUpperBound(1);
    //params->setStepsizeAccuracy(0.0005);
    //params->setVerbosity("DETAIL");
    //params->enableLog();
    //optimizer.init(initial_x, f, params);

    //MomentumOptimizer<RosenbrockFunctor> optimizer;
    //MomentumParams* params = new MomentumParams;
    //params->setVerbosity("DETAIL");
    //params->enableLog();
    //optimizer.init(initial_x, f, params);

    //NesterovMomentumOptimizer<RosenbrockFunctor> optimizer;
    //NesterovMomentumParams* params = new NesterovMomentumParams;
    //params->setVerbosity("DETAIL");
    //params->enableLog();
    //optimizer.init(initial_x, f, params);

    //AdagradOptimizer<RosenbrockFunctor> optimizer;
    //AdagradParams* params = new AdagradParams;
    //params->enableLog();
    //params->setVerbosity("DETAIL");
    //optimizer.init(initial_x, f, params);

    //RMSPropOptimizer<RosenbrockFunctor> optimizer;
    //RMSPropParams* params = new RMSPropParams;
    //params->setVerbosity("DETAIL");
    //params->enableLog();
    //optimizer.init(initial_x, f, params);

    //AdaDeltaOptimizer<RosenbrockFunctor> optimizer;
    //AdaDeltaParams* params = new AdaDeltaParams;
    //params->enableLog();
    //params->setSGDTimes(1100);
    //params->setVerbosity("DETAIL");
    //optimizer.init(initial_x, f, params);

    //AdamOptimizer<RosenbrockFunctor> optimizer;
    //AdamParams* params = new AdamParams;
    //params->enableLog();
    //optimizer.init(initial_x, f, params);

    //HookeJeevesOptimizer<RosenbrockFunctor> optimizer;
    //HookeJeevesParams* params = new HookeJeevesParams;
    //params->setVerbosity("DETAIL");
    //params->enableLog();
    //optimizer.init(initial_x, f, params);

    //NewtonOptimizer<RosenbrockFunctor, RosenbrockHessianFunctor> optimizer;
    //NewtonParams* params = new NewtonParams;
    //params->enableLog();
    //params->setVerbosity("DETAIL");
    //optimizer.init(initial_x, f, params);

    //LevenbergMarquardtOptimizer<RosenbrockFunctor, RosenbrockHessianFunctor> optimizer;
    //LevenbergMarquardtParams* params = new LevenbergMarquardtParams;
    //params->setVerbosity("DETAIL");
    //params->enableLog();
    //optimizer.init(initial_x, f, params);

    //DFP_Optimizer<RosenbrockFunctor> optimizer;
    //DFP_Params* params = new DFP_Params;
    //params->setStepsizeAccuracy(0.002);
    //params->setStepsizeIncreaseFactor(2);
    //params->setStepsizeDecreaseFactor(0.3);
    //params->enableLog();
    //params->setVerbosity("DETAIL");
    //optimizer.init(initial_x, f, params);

    BFGS_Optimizer<RosenbrockFunctor> optimizer;
    BFGS_Params* params = new BFGS_Params;
    params->enableLog();
    params->setVerbosity("DETAIL");
    params->setStepsizeMethod("GOLDENSECTION");
    optimizer.init(initial_x, f, params);

    optimizer.optimize();

    return 0;
}

