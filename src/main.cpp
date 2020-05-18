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

struct PhiFunctor: public NLOP::Functor<double, 1>
{
    double a = 1;
    double b = 100;

    Eigen::Matrix<double, 1, 2> J;
    Eigen::Matrix<double, 2, 1> P;

    void setJ(const Eigen::Matrix<double, 1, 2>& value)
    {
        J = value;
    }

    void setP(const Eigen::Matrix<double, 2, 1>& point)
    {
        P = point;
    }

    double operator() (double lambda)
    {
        //double lambda = x[0];
        double J_x1 = J[0];
        double J_x2 = J[1];
        double x1 = P[0];
        double x2 = P[1];

        return (a - (x1 - lambda*J_x1))*(a - (x1 - lambda*J_x1))
                + b*((x2 - lambda*J_x2) - (x1 - lambda*J_x1)*(x1 - lambda*J_x1))
                *((x2 - lambda*J_x2) - (x1 - lambda*J_x1)*(x1 - lambda*J_x1));
    }
};

using namespace NLOP;

int main()
{

    /*
    // Steepest Descent Optimization Test
    RosenbrockFunctor* f = new RosenbrockFunctor;
    SteepestDescentParams* params= new SteepestDescentParams;
    Eigen::Vector2d initial_x(0, 0);
    SteepestDescentOptimizer<double, 2, RosenbrockFunctor, PhiFunctor> optimizer;
    optimizer.init(initial_x, f, params);
    optimizer.optimize();
    */

    // Conjuate Gradient Optimization Test
    RosenbrockFunctor* f = new RosenbrockFunctor;
    ConjuateGradientParams* params = new ConjuateGradientParams;
    Eigen::Vector2d Initial_x(0, 0);
    ConjuateGradientOptimizer<double, 2, RosenbrockFunctor, PhiFunctor> optimizer;
    params->setApproximateMethod("PR");
    optimizer.init(Initial_x, f, params);
    optimizer.optimize();

    return 0;
}

