cmake_minimum_required(VERSION 2.8)

project(ConvexOptimization)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "-std=c++11 -O0")

find_package(Eigen3 REQUIRED)

include_directories(
    "include"
    ${Eigen3_INCLUDE_DIR}
)

set(HEADERS
    # Utils
    include/Utils/Matrix.hpp
    include/Utils/Types.hpp
    include/Utils/Utils.hpp
    include/Utils/Functor.hpp

    # Optimizer Base
    include/OptimizerBase/OptimizerBase.hpp
    include/OptimizerBase/LineSearchOptimizer.hpp
    include/OptimizerBase/TrustRegionOptimizer.hpp

    #Optimizer Params
    include/OptimizerParams/OptimizerParamsBase.hpp
    include/OptimizerParams/LineSearchParams.hpp
    include/OptimizerParams/TrustRegionParams.hpp
    include/OptimizerParams/SteepestDescentParams.hpp
    include/OptimizerParams/ConjuateGradientParams.hpp

    # Line Search Methods
    include/LineSearchMethods/CyclicCoordinateOptimizer.hpp
    include/LineSearchMethods/SteepestDescentOptimizer.hpp
    include/LineSearchMethods/NewtonOptimizer.hpp
    include/LineSearchMethods/ConjuateGradientOptimizer.hpp

    # Trust Region Methods
    include/TrustRegionMethods/LevenbergMarquardtOptimizer.hpp

    # One Dimensional Search Methods
    include/OneDimensionalSearch/OneDimensionalSearchMethods.hpp

    include/OneDimensionalSearch/GoldenSectionMethod.hpp
    include/OneDimensionalSearch/FibonacciMethod.hpp
    include/OneDimensionalSearch/BisectionMethod.hpp
    include/OneDimensionalSearch/DichotomousMethod.hpp
    include/OneDimensionalSearch/NewtonMethod.hpp
    )

# Examples
#if(CXX_BUILD_EXAMPLES)
#  add_executable(example_robot1 "examples/Robot1/main.cpp")
#endif()

add_executable(${PROJECT_NAME}
    "src/main.cpp")
