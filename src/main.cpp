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

//#include <ConstrainedOptimizer/ConstrainedOptimizerBase.hpp>

using namespace std;

struct RosenbrockFunctor: public NLOP::Functor<double, 2>
{
    //template <typename T1, typename T2, typename T3>
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
    //params->setStepsizeMethod("ARMIJO");
    //optimizer.init(initial_x, f, params);

    //ConjuateGradientOptimizer<RosenbrockFunctor> optimizer;
    //ConjuateGradientParams* params = new ConjuateGradientParams;
    //optimizer.init(initial_x, f, params);

    //MomentumOptimizer<RosenbrockFunctor> optimizer;
    //MomentumParams* params = new MomentumParams;
    //optimizer.init(initial_x, f, params);

    //NesterovMomentumOptimizer<RosenbrockFunctor> optimizer;
    //NesterovMomentumParams* params = new NesterovMomentumParams;
    //optimizer.init(initial_x, f, params);

    //AdagradOptimizer<RosenbrockFunctor> optimizer;
    //AdagradParams* params = new AdagradParams;
    //optimizer.init(initial_x, f, params);

    //RMSPropOptimizer<RosenbrockFunctor> optimizer;
    //RMSPropParams* params = new RMSPropParams;
    //optimizer.init(initial_x, f, params);

    //AdaDeltaOptimizer<RosenbrockFunctor> optimizer;
    //AdaDeltaParams* params = new AdaDeltaParams;
    //optimizer.init(initial_x, f, params);

    //AdamOptimizer<RosenbrockFunctor> optimizer;
    //AdamParams* params = new AdamParams;
    //optimizer.init(initial_x, f, params);

    //HookeJeevesOptimizer<RosenbrockFunctor> optimizer;
    //HookeJeevesParams* params = new HookeJeevesParams;
    //params->setVerbosity("DETAIL");
    //optimizer.init(initial_x, f, params);

    //NewtonOptimizer<RosenbrockFunctor, RosenbrockHessianFunctor> optimizer;
    //NewtonParams* params = new NewtonParams;
    //optimizer.init(initial_x, f, params);

    LevenbergMarquardtOptimizer<RosenbrockFunctor, RosenbrockHessianFunctor> optimizer;
    LevenbergMarquardtParams* params = new LevenbergMarquardtParams;
    optimizer.init(initial_x, f, params);

    //DFP_Optimizer<RosenbrockFunctor> optimizer;
    //DFP_Params* params = new DFP_Params;
    //optimizer.init(initial_x, f, params);

    //BFGS_Optimizer<RosenbrockFunctor> optimizer;
    //BFGS_Params* params = new BFGS_Params;
    //optimizer.init(initial_x, f, params);

    optimizer.optimize();

    return 0;
}

