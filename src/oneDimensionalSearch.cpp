#include <iostream>
#include <oneDimensionalSearch.h>

using namespace std;

/// Golden Section Search Method
/// interval: initial searching interval
/// t: reduction ratio (default: 0.618)
/// epsilon: target width of interval
double goldenSectionSearch(std::pair<double, double> interval, double epsilon, double t = 0.618)
{
    // initialization
    double alpha = interval.first;
    double beta = interval.second;
    double lambda = alpha + (1 - t) * (beta - alpha);
    double mu = alpha + t * (beta - alpha);
    int iterative_count = 0;

    // iteratively narrow the searching interval to find the optimal stepsize
    while(true)
    {
        cout << "interative times: " << iterative_count << "  "
             << "[alpha, lambda, mu, beta]: "
             << "[" << alpha << ", " << lambda << ", "
             << mu << ", " << beta << "]" << endl;

        // stoping condition: (alpha - beta) < epsilon
        if (beta - alpha < epsilon)
        {
            cout << "Finished! Optimal lambda: "
                 << (alpha + beta)/2 << endl;
            return (alpha + beta)/2;
        }
        iterative_count++;

        if (phi(lambda) > phi(mu))
        //if (f(lambda) > f(mu))
        {
            // cut the interval [alpha, lambda)
            alpha = lambda;
            lambda = mu;
            mu = alpha + t * (beta - alpha);
        }
        else
        {
            // cut the interval (mu, beta]
            beta = mu;
            mu = lambda;
            lambda = alpha + (1 - t) * (beta - alpha);
        }
    }
}



/// Fibonacci Search Method
/// interval: initial searching interval
/// n: the index of used max fibonacci number, determined on required accuracy
/// epsilon: target width of interval
double fibonacciSearch(std::pair<double, double> interval, double epsilon, int n)
{
    // initialization
    double alpha = interval.first;
    double beta = interval.second;
    double t = fibonacci(n-1)/fibonacci(n);
    double lambda = alpha + (1 - t) * (beta - alpha);
    double mu = alpha + t * (beta - alpha);
    int iterative_count = 0;

    // iteratively narrow the searching interval with the varying ratio to find the optimal stepsize
    while (true)
    {
        cout << "interative times: " << iterative_count << "  "
             << "[alpha, lambda, mu, beta]: "
             << "[" << alpha << ", " << lambda << ", "
             << mu << ", " << beta << "], t: " << t << endl;

        // stoping condition: (alpha - beta) < epsilon
        if (beta - alpha < epsilon)
        {
            cout << "Finished! Optimal lambda: "
                 << (beta + alpha)/2 << endl;
            return (beta + alpha)/2;
        }
        ++iterative_count;

        //if (phi(lambda) > phi(mu))
        if (f(lambda) > f(mu))
        {
            alpha = lambda;
            lambda = mu;
            mu = alpha + t * (beta - alpha);
        }
        else
        {
            beta = mu;
            mu = lambda;
            lambda = alpha + (1 - t) * (beta - alpha);
        }

        if (--n == 0)
        {
            cerr << "n is too small, please choose a bigger n and try again" << endl;
            return -1;
        }
        // update reduction ratio;
        t = fibonacci(n-1)/fibonacci(n);
    }
}

/// Dichotomous Method for General Continues Functions
/// interval: initial searching interval
/// epsilon: 1/2 target width of interval
double dichotomousSearch(std::pair<double, double> interval, double epsilon)
{
    // initialization
    double alpha = interval.first;
    double beta = interval.second;
    double lambda = (alpha + beta)/2 - epsilon;
    double mu = (alpha + beta)/2 + epsilon;
    int iterative_count = 0;

    int max_iteration = 12;

    // iteratively cut a half of the interval to find the optimal lambda
    while (true)
    {
        cout << "interative times: " << iterative_count << "  "
             << "[alpha, lambda, mu, beta]: "
             << "[" << alpha << ", " << lambda << ", "
             << mu << ", " << beta << "]" << endl;

        // stoping condition: (alpha - beta) < 2*epsilon or reach max iteration times
        if (beta - alpha < 2*epsilon || iterative_count == max_iteration)
        {
            cout << "Finished! Optimal lambda: "
                 << lambda << endl;
            return lambda;
        }
        ++iterative_count;

        //if (phi(mu) > phi(lambda))
        if (f(mu) > f(lambda))
        {
            // cut the right half of the interval
            beta = mu;
            lambda = (alpha + beta)/2 - epsilon;
            mu = (alpha + beta)/2 + epsilon;
        }

        else {
            // cut the left half of the interval
            alpha = lambda;
            lambda = (alpha + beta)/2 - epsilon;
            mu = (alpha + beta)/2 + epsilon;
        }
    }
}

/// Dichotomous Method for Derivable Functions
/// interval: initial searching interval
/// epsilon: target width of interval
double bisectionSearch(std::pair<double, double> interval, double epsilon)
{
    // initialization
    double alpha = interval.first;
    double beta = interval.second;
    double lambda = (alpha + beta)/2;
    int iterative_count = 0;

    // iteratively cut a half of the interval to find the optimal lambda
    while (true)
    {
        cout << "interative times: " << iterative_count << "  "
             << "[alpha, lambda, beta]: "
             << "[" << alpha << ", " << lambda << ", "
             << beta << "]" << endl;

        // stoping condition: (alpha - beta) < epsilon || phi'(lambda) == 0
        //if (beta - alpha < epsilon || phiDerivative(lambda) == 0)
        if (beta - alpha < epsilon || fDerivative(lambda) == 0)
        {
            cout << "Finished! Optimal lambda: "
                 << lambda << endl;
            return lambda;
        }
        ++iterative_count;

        //if (phiDerivative(lambda) > 0)
        if (fDerivative(lambda) > 0)
        {
            // cut the right half of the interval
            beta = lambda;
            lambda = (alpha + beta)/2;
        }

        else
        {
            // cut the left half of the interval
            alpha = lambda;
            lambda = (alpha + beta)/2;
        }
    }
}

double goldenSectionSearchRosenBrock(std::pair<double, double> interval, double epsilon,
                                     Eigen::Vector2d recent_point, Eigen::Vector2d Jacobian,
                                     double a = 1, double b = 100, double t = 0.618)
{
    // initialization
    double alpha = interval.first;
    double beta = interval.second;
    double lambda = alpha + (1 - t) * (beta - alpha);
    double mu = alpha + t * (beta - alpha);
    int iterative_count = 0;

    // iteratively narrow the searching interval to find the optimal stepsize
    while(true)
    {
        cout << "interative times: " << iterative_count << "  "
             << "[alpha, lambda, mu, beta]: "
             << "[" << alpha << ", " << lambda << ", "
             << mu << ", " << beta << "]" << endl;

        // stoping condition: (alpha - beta) < epsilon
        if (beta - alpha < epsilon)
        {
            cout << "Finished! Optimal lambda: "
                 << (alpha + beta)/2 << endl;
            return (alpha + beta)/2;
        }
        iterative_count++;

        if (phiRosenBrock(recent_point, Jacobian, lambda, a, b) > phiRosenBrock(recent_point, Jacobian, mu, a, b))
        {
            // cut the interval [alpha, lambda)
            alpha = lambda;
            lambda = mu;
            mu = alpha + t * (beta - alpha);
        }
        else
        {
            // cut the interval (mu, beta]
            beta = mu;
            mu = lambda;
            lambda = alpha + (1 - t) * (beta - alpha);
        }
    }
}


