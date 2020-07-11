NLOP - Non-Linear Optimization Library
===================================================

# Introduction

NLOP is an open-source C++ template library for non-linear optimization.


# Motivation

I developed this project mainly because of my interest in convex optimization.
I hope this will be useful for other convex optimization researchers for learning and teaching.

# Features

- Thanks to the template programming, users can customize the data type to trade-off between speed and precision
- Autodiff function is provided, users don't have to give Jacobian formula manually
- Line search methods can be combined with step search method flexibly
- Users can enable log to have the optimization process information saved in an txt file, so that they can be used for visualization (a python script example is provided)

# Optimization method list

### Linesearch methods
- Hooke & Jeeves's method
- steepest descent method
- conjuate gradient method
- Momentum
- Nesterov Momentum
- Adagrad
- RMSPro
- AdaDelta
- Adam

### Newton's and quasi-Newton methods
- Newton's method
- DFP
- BFGS

### Trust region methods
- Levenberg-Marquardt method

Stepsize searching method list
------------------------------

## Accurate searching methods
- Golden section method
- Fibonacci method
- Dichotomous method

## Inaccurate searching methods
- Armijo method
- Goldstein method
- Wolfe-Powell method

# Requirements

- [cmake](http://www.cmake.org/)
- [Eigen3](http://eigen.tuxfamily.org)

# Contacts
--------
Please use Github issues to report bugs. If you have ang question or advice, please feel free to contact [Shuang Guo](mailto:guoshuangSLAM@outlook.com).
