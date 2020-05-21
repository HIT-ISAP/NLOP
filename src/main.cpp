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
#include <LineSearchMethods/RosenbrockOptimizer.hpp>

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

struct HookeTestFunctor: public NLOP::Functor<double, 2>
{
    double x1, x2;
    double operator ()(const InputType& x)
    {
        x1 = x[0];
        x2 = x[1];
        double result = pow(x1 - 2, 4) + pow(x1 - 2*x2, 2);
        return result;
    }

    void operator ()(const ActiveInput& x, ActiveValue* v) const
    {
        (*v)[0] = pow(x1 - 2, 4) + pow(x1 - 2*x2, 2);
    }

    void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const
    {
        (*v)[0] = pow(x[0] - 2, 4) + pow(x[0] - 2*x[1], 2);
        std::cout << (*v)[0] << std::endl;
    }
};

using namespace NLOP;

int main()
{
    RosenbrockFunctor* f = new RosenbrockFunctor;
    Eigen::Vector2d initial_x(0, 3);

    //auto ss = new GoldSectionMethod<RosenbrockFunctor>;
    //auto ss = new FibonacciMethod<RosenbrockFunctor>;
    //auto ss = new DichotomousMethod<RosenbrockFunctor>;

    //auto ss = new GoldsteinMethod<RosenbrockFunctor>;
    //auto ss = new ArmijoMethod<RosenbrockFunctor>;
    //auto ss = new WolfePowellMethod<RosenbrockFunctor>;

    //SteepestDescentOptimizer<RosenbrockFunctor> optimizer;
    //SteepestDescentParams* params= new SteepestDescentParams;

    //ConjuateGradientOptimizer<RosenbrockFunctor> optimizer;
    //ConjuateGradientParams* params = new ConjuateGradientParams;

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

    HookeTestFunctor* f_t = new HookeTestFunctor;
    HookeJeevesOptimizer<RosenbrockFunctor> optimizer;
    HookeJeevesParams* params = new HookeJeevesParams;
    optimizer.init(initial_x, f, params);

    //optimizer.init(initial_x, f, params, ss);
    optimizer.optimize();

    return 0;
}

