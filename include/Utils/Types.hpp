#ifndef NLOP_TYPES_HPP
#define NLOP_TYPES_HPP

#include <Utils/Matrix.hpp>

namespace NLOP {

/// @class NLOP::SquareMatrix
/// @brief Template type representing a square matrix
/// @param T The numeric scalar type
/// @param N The dimension of the square matrix
template<typename T, int N>
using SquareMatrix = Matrix<T, N, N>;

/// @class NLOP::Jacobian
/// @brief Template type of Jacobian of 1 w.r.t N
/// @param VariableType Variable x with a dimension of N
// template <typename VariableType>
// using Jacobian = Vector<typename VariableType::Scalar, typename VariableType::RowsAtCompileTime>;

/// @class NLOP::Hessian
/// @brief Template type of Hessian of dimension N
/// @param VariableType Variable x with a dimension of
// template <typename VariableType>
// using Hessian = SquareMatrix<typename VariableType::Scalar, typename VariableType::RowsAtCompileTime>;

}

#endif
