#include <iostream>
#include <Utils/Functor.hpp>
#include <LineSearchMethods/SteepestDescentOptimizer.hpp>
#include <LineSearchMethods/ConjuateGradientOptimizer.hpp>

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


using namespace NLOP;

int main()
{
    RosenbrockFunctor* f = new RosenbrockFunctor;
    Eigen::Vector2d initial_x(0, 0);

    //auto ss = new GoldSectionMethod<RosenbrockFunctor>;
    //auto ss = new FibonacciMethod<RosenbrockFunctor>;
    //auto ss = new DichotomousMethod<RosenbrockFunctor>;

    auto ss = new GoldsteinMethod<RosenbrockFunctor>;

    //SteepestDescentOptimizer<RosenbrockFunctor> optimizer;
    //SteepestDescentParams* params= new SteepestDescentParams;

    ConjuateGradientOptimizer<RosenbrockFunctor> optimizer;
    ConjuateGradientParams* params = new ConjuateGradientParams;

    optimizer.init(initial_x, f, params, ss);
    optimizer.optimize();



    return 0;
}

