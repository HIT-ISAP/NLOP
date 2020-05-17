#ifndef NLOP_FUNCTOR_HPP
#define NLOP_FUNCTOR_HPP

#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/AutoDiff>

namespace NLOP {

/// @class NLOP::Functor
/// @brief Abstract template type for function expression
/// @param T The numeric scalar type
/// @param N the dimension of input variable x
template<typename T, int N>
class Functor
{
public:
    using InputType = Eigen::Matrix<T, N, 1>;
    using ValueType = Eigen::Matrix<T, 1, 1>;
    using JacobianType = Eigen::Matrix<T, 1, N>;

    enum{
        InputsAtCompileTime = InputType::RowsAtCompileTime,
        ValuesAtCompileTime = ValueType::RowsAtCompileTime
    };

    using DerivativeType = Eigen::Matrix<T,InputsAtCompileTime,1>;
    using ActiveScalar = Eigen::AutoDiffScalar<DerivativeType>;

    using ActiveInput = Eigen::Matrix<ActiveScalar, InputsAtCompileTime, 1>;
    using ActiveValue = Eigen::Matrix<ActiveScalar, ValuesAtCompileTime, 1>;

    int m_inputs = InputsAtCompileTime;
    int m_values = ValuesAtCompileTime;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime)
    {
        y = new ValueType;
        jacobian = new JacobianType;
    }
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values)
    {
        y = new ValueType;
        jacobian = new JacobianType;
    }

    virtual ~Functor()
    {
        delete y;
        delete jacobian;
    }

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

    T getY()
    {
        return (*y)(0,0);
    }

    JacobianType getJacobian()
    {
        return (*jacobian);
    }

    InputType x;
    ValueType* y;
    JacobianType* jacobian;

    // you should define that in the subclass :
    // void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
    // void operator() (const ActiveInput& x, ActiveValue* v) const;
};
}

#endif
