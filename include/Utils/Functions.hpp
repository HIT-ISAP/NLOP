#ifndef NLOP_FUNCTIONS_HPP
#define NLOP_FUNCTIONS_HPP

#include <Utils/Types.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/AutoDiff>

namespace NLOP {

/// @class NLOP::Function
/// @brief Template type representing a non-linear function
/// @param T The numeric scalar type
/// @param N The dimension of variable x
template<typename VariableType>
class Functor
{
public:

    /// Variable type
    using Variable = VariableType;

    /// Variable dimension
    //using N = typename Variable::RowsAtCompileTime;
    //using M = typename OutputType::RowsAtCompileTime;

    /// The numeric scalar type of value
    using T = typename VariableType::Scalar;


    using InputType = Eigen::Matrix<T, 3, 1>;
    using ValueType = Eigen::Matrix<T, 1, 1>;
    using JacobianType = Eigen::Matrix<T, 1, 3>;

    enum {
        InputsAtCompileTime = InputType::RowsAtCompileTime,
        ValuesAtCompileTime = ValueType::RowsAtCompileTime
      };

    typedef Eigen::Matrix<T,InputsAtCompileTime,1> DerivativeType;
    typedef Eigen::AutoDiffScalar<DerivativeType> ActiveScalar;


    typedef Eigen::Matrix<ActiveScalar, InputsAtCompileTime, 1> ActiveInput;
    typedef Eigen::Matrix<ActiveScalar, ValuesAtCompileTime, 1> ActiveValue;

    int m_inputs = InputsAtCompileTime, m_values = ValuesAtCompileTime;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; } // number of degree of freedom (= 2*nb_vertices)
    int values() const { return m_values; } // number of energy terms (= nb_vertices + nb_edges)


    void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const
    {
        (*v)[0] = log(x[0]) + log(x[1]) + log(x[2]);
    }

    void operator ()(const ActiveInput& x, ActiveValue* v) const
    {
        (*v)[0] = log(x[0]) + log(x[1]) + log(x[2]);
    }

    T operator() (const InputType& x) const
    {
        operator ()(x, value);
        return value(0, 0);
    }

    /*
    /// @brief Compute the value of f(x) with given x
    T f(const Variable& x) const
    {
        ValueType temp_value;
        operator ()(x, temp_value);
        return temp_value(0, 0);
    }

    /// @brief Compute the value of f(x) with current x

    T f()
    {
        operator ()(x, value);
        return value(0, 0);
    }
    */

    /// @brief Set the value of x
    void setX(const Variable& x)
    {
        this->x = x;
    }

    /// @brief Get the value of current x
    T getX()
    {
        return x;
    }

    /// @brief update Jacobian matrix of f(x)
    JacobianType updateJacobian()
    {
        //adjac(x, &value, &Jacobian);
        return this->Jacobian;
    }

    /// @brief get J(x) with current x
    JacobianType getJacobian()
    {
        return this->Jacobian;
    }

    /*
    /// @brief update Hessian matrix of f(x)
    void updateHessian()
    {
        /// @todo
    }

    /// @brief get H(x) with current x
    Matrix getHessian()
    {
        return this->Hessian;
    }
    */

    //Functor() {}
    //Functor(Variable x): x(x) {isInit = true;}
    //~Fucntion() {}

//private:
    ValueType value; // Output Value
    //ValueType ratio; // The ratio using for scalarization

    Variable x; // Variable x

    JacobianType Jacobian; // Jacobian Matrix


    /// H(x)
    // Matrix Hessian;

    /// @brief Function f(x), compute f(x) with given x
    //T f(const Variable& x);
    /// TODO
};
}



#endif
