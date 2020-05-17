#ifndef NEWTONMETHOD_HPP
#define NEWTONMETHOD_HPP

#include <OneDimensionalSearch/OneDimensionalSearchMethods.hpp>

namespace NLOP {

/// @class NLOP::NewtonMethod
/// @brief Newton method for second order differentiable function
/// @param T The numeric scalar type
template <typename T>
class NewtonMethod: public OneDimSearch
{
public:
    T search() override
    {
        /// TODO
    }

};
}
