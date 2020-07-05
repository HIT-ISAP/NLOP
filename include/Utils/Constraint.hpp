#ifndef CONSTRAINT_HPP
#define CONSTRAINT_HPP

#include <eigen3/Eigen/Core>

namespace NLOP {
/// @class NLOP::EqualityConstraint
/// @brief Abstract template type for linear constraint expression: Ax = b
/// @param T The numeric scalar type
/// @param M The dimension of the rhs vector b
/// @param N The dimension of input variable x
template<typename T, int M, int N>
class LinearConstraint
{
public:
    using InputType = Eigen::Matrix<T, N, 1>;
    using ResidualType = Eigen::Matrix<T, M, 1>;
    using CoeffMatrixType = Eigen::Matrix<T, M, N>;
    using Scalar = T;

    enum{
        InputsAtCompileTime = InputType::RowsAtCompileTime,
        ResidualsAtCompileTime = ResidualType::RowsAtCompileTime
    };

    int m_inputs = InputsAtCompileTime;
    int m_values = ResidualsAtCompileTime;

    LinearConstraint(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    virtual ~LinearConstraint() {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

    void setX(const InputType& x_value)
    {
        x = x_value;
    }

    InputType getX()
    {
        return x;
    }

    InputType x;
    ResidualType b;
    CoeffMatrixType A;
};
}

#endif
