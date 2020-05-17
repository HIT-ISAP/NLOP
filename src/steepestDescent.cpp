#include <steepestDescent.h>
#include <fstream>

Eigen::Vector2d steepestDescent(Eigen::Vector2d& start_point, double epsilon)
{
    // ofstream output("iteration.txt");
    Eigen::Vector2d Jacobian; // inverse direction of gradient
    Eigen::Vector2d recent_point = start_point;
    double lambda = 0; // stepsize

    std::pair<double, double> initial_interval(0, 1);

    int k = 1; // kth iterative point
    double recent_value = 0;

    // output << recent_point[0] << " " << recent_point[1] << " " << recent_value << "\n";

    while (true)
    {
        // compute gradient
        Jacobian = rosenBrockDerivative(recent_point, 1, 100);

        if (Jacobian.norm() < epsilon)
        {
            cout << "gradient: " << Jacobian.norm()
                 << ", Finished! " << endl;
            return recent_point;
        }

        else
        {
            k++;
            // cout << "gradient norm: " << Jacobian.norm() << endl;

            // search optimal stepsize lambda using golden section method
            lambda = goldenSectionSearchRosenBrock(initial_interval, 0.001, recent_point, Jacobian, 1, 100, 0.618);
            // update recent point and value
            recent_point += lambda * (-Jacobian);
            recent_value = rosenBrock(recent_point, 1, 100);

            // output << recent_point[0] << " " << recent_point[1] << " " << recent_value << "\n";

            cout << k << " iteration: " << recent_point << ", "
                 << "recent value: " << recent_value << endl;

        }
    }
}
