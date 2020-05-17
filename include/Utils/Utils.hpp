#ifndef NLOP_UTILS_H
#define NLOP_UTILS_H

#include <Utils/Matrix.hpp>
#include <iostream>

namespace NLOP
{
/// @brief Recursively calculate Fibonacci number
/// @param n Fibonacci array index
double fibonacci(int n)
{
    if (n < 1)
    {
        std::cerr << "Invalid index n!" << std::endl;
        return -1;
    }
    if (n == 1 || n == 2)
        return 1;
    else
        return fibonacci(n-1) + fibonacci(n-2);
}

}

#endif
