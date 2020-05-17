#include <utils.h>

using namespace std;

/// objective function related to stepsize lambda
double phi(double lambda)
{
    return lambda*lambda + 4*lambda + 5;
}
/// deriviative for function phi
double phiDerivative(double lambda)
{
    return 2*lambda + 4;
}

/// function f(x) in question 13 of homework
double f(double x)
{
    return x*x/2 - x;
}

double fDerivative(double x)
{
    return x-1;
}

/// recursively calculate Fibonacci number for index n
double fibonacci(int n)
{
    if (n < 1)
    {
        cerr << "Invalid index n!" << endl;
        return -1;
    }
    if (n == 1 || n == 2)
        return 1;
    else
        return fibonacci(n-1) + fibonacci(n-2);
}

/// Rosenbrock function: (a - x1)^2 + b(x2 - x1^2)^2
/// params a, b is variant
double rosenBrock(Eigen::Vector2d x, double a = 1, double b = 100)
{
    double x1 = x[0];
    double x2 = x[1];
    return (a - x1)*(a - x1) + b*(x2 - x1*x1)*(x2 - x1*x1);
}

/// Jacobian for Rosenbrock function
Eigen::Vector2d rosenBrockDerivative(Eigen::Vector2d x, double a = 1, double b = 100)
{
    Eigen::Vector2d Jacobian;
    double x1 = x[0];
    double x2 = x[1];
    Jacobian[0] = 4*b*x1*x1*x1 - 4*b*x1*x2 + 2*x1 - 2*a;
    Jacobian[1] = 2*b*(x2 - x1*x1);
    return Jacobian;
}

/// function phi for Rosenbrock function
double phiRosenBrock(Eigen::Vector2d x, Eigen::Vector2d Jacobian,
                     double lambda, double a = 1, double b = 100)
{
    double x1 = x[0];
    double x2 = x[1];
    double J_x1 = Jacobian[0];
    double J_x2 = Jacobian[1];

    return (a - (x1 - lambda*J_x1))*(a - (x1 - lambda*J_x1))
            + b*((x2 - lambda*J_x2) - (x1 - lambda*J_x1)*(x1 - lambda*J_x1))
            *((x2 - lambda*J_x2) - (x1 - lambda*J_x1)*(x1 - lambda*J_x1));
}
