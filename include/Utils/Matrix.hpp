#ifndef NLOP_MATRIX_HPP
#define NLOP_MATRIX_HPP

#include <cmath>

#include <Eigen/Dense>

namespace NLOP {

/// @class NLOP::Matrix
/// @brief Template class for matrices
/// @param T The numeric scalar type
/// @param rows The number of rows
/// @param cols The number of columns
template<typename T, int rows, int cols>
using Matrix = Eigen::Matrix<T, rows, cols>;

/// @class NLOP::Vector
/// @brief Template class for vectors
/// @param T The numberic scalar type
/// @param N The vector dimension
template<typename T, int N>
class Vector: public Matrix<T, N, 1>
{
public:
    typedef Matrix<T, N, 1> Base;
    using typename Base::Scalar;
    using Base::RowsAtCompileTime;
    using Base::ColsAtCompileTime;
    using Base::SizeAtCompileTime;

    /// @brief Constructor
    Vector(void): Matrix<T, N, 1>() {}

    /// @brief Copy constructor
    template<typename OtherDerived>
    Vector(const Eigen::MatrixBase<OtherDerived>& other)
        : Matrix<T, N, 1>(other) {}

    /// @brief Copy assignment constructor
    template<typename OtherDerived>
    Vector& operator =(const Eigen::MatrixBase <OtherDerived>& other)
    {
        this->Base::operator =(other);
        return *this;
    }
};


}

#endif
