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

/// the definition of Rosenbrock function, for your own obejective function, please imitate this example
/// note that these 3 forms of operator() are all necessary to ensure all the optimization and stepsize searching methods to run smoothly
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

/// the auto-calculation for Hessian matrix is not available at the moment
/// so if Newton's method or L-M method is used, a Hessian matrix functor is needed
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
    RosenbrockFunctor* f = new RosenbrockFunctor; /// initialize a pointer to the objective function,
                                                  /// you don't have to delete this pointer, this has been done in the deconstructor
    Eigen::Vector2d initial_x(0, 0);              /// initial x value to launch the optimization

    SteepestDescentOptimizer<RosenbrockFunctor> optimizer;      /// initialize a steepest descent optimizer, other optimizer are given below,
                                                                /// you can replace this to try them
    SteepestDescentParams* params= new SteepestDescentParams;   /// initialize default parameters for steepest descent optimizer,
                                                                /// you can modify some of them by using the setters, examples are given below
                                                                /// for the details of all the interface, see their param class headers

    //params->setVerbosity("DETAIL");             // choose the mode for information printing
    //params->setMaxItertaions(1000);             // set max iteration times for optimization
    //params->setMinGradient(0.1);                // set gradient threshold to stop the iterations
    //params->setStepsizeMethod("WOLFEPOWELL");   // set stepsize searching method, 6 types are available now (3 accurate types and 3 inaccurate types)
    //params->setStepsizeAccuracy(0.001);         // set stepsize accuracy (only for accurate searching methods)
    //params->setStepsizeIncreaseFactor(3);       // set stepsize increase factor (only for inaccurate searching methods)
    //params->setStepsizeDecreaseFactor(0.3);     // set stepsize decrease factor (only for inaccurate searching methods)
    //params->setStepsizeUpperBound(0.01);        // set stepsize upper bound
    //params->enableLog();                        // enable to write optimization process information to a txt file, for plotting purpose


    // other available optimizers
    //ConjuateGradientOptimizer<RosenbrockFunctor> optimizer;
    //ConjuateGradientParams* params = new ConjuateGradientParams;

    //MomentumOptimizer<RosenbrockFunctor> optimizer;
    //MomentumParams* params = new MomentumParams;

    //NesterovMomentumOptimizer<RosenbrockFunctor> optimizer;
    //NesterovMomentumParams* params = new NesterovMomentumParams;

    //AdagradOptimizer<RosenbrockFunctor> optimizer;
    //AdagradParams* params = new AdagradParams;

    //RMSPropOptimizer<RosenbrockFunctor> optimizer;
    //RMSPropParams* params = new RMSPropParams;

    //AdaDeltaOptimizer<RosenbrockFunctor> optimizer;
    //AdaDeltaParams* params = new AdaDeltaParams;

    //AdamOptimizer<RosenbrockFunctor> optimizer;
    //AdamParams* params = new AdamParams;

    //HookeJeevesOptimizer<RosenbrockFunctor> optimizer;
    //HookeJeevesParams* params = new HookeJeevesParams;

    //NewtonOptimizer<RosenbrockFunctor, RosenbrockHessianFunctor> optimizer;
    //NewtonParams* params = new NewtonParams;

    //LevenbergMarquardtOptimizer<RosenbrockFunctor, RosenbrockHessianFunctor> optimizer;
    //LevenbergMarquardtParams* params = new LevenbergMarquardtParams;

    //DFP_Optimizer<RosenbrockFunctor> optimizer;
    //DFP_Params* params = new DFP_Params;

    //BFGS_Optimizer<RosenbrockFunctor> optimizer;
    //BFGS_Params* params = new BFGS_Params;

    optimizer.init(initial_x, f, params);   /// initialize the optimizer with initial_value, objective function and given params
    optimizer.optimize();                   /// start the optimization process

    delete params;

    return 0;
}

