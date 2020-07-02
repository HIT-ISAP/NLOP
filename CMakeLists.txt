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
    include/Utils/Hessian.hpp

    #Optimizer Params Base
    include/OptimizerParams/OptimizerParamsBase.hpp
    include/OptimizerParams/LineSearchParams.hpp
    include/OptimizerParams/TrustRegionParams.hpp

    #Line search method params
    include/OptimizerParams/SteepestDescentParams.hpp
    include/OptimizerParams/ConjuateGradientParams.hpp
    include/OptimizerParams/MomentumParams.hpp
    include/OptimizerParams/NesterovMomentumParams.hpp
    include/OptimizerParams/AdagradParams.hpp
    include/OptimizerParams/RMSPropParams.hpp
    include/OptimizerParams/AdaDeltaParams.hpp
    include/OptimizerParams/AdamParams.hpp
    include/OptimizerParams/HookeJeevesParams.hpp
    include/OptimizerParams/RosenbrockParams.hpp

    include/OptimizerParams/NewtonParams.hpp
    include/OptimizerParams/LevenbergMarquardtParams.hpp
    include/OptimizerParams/DFP_Params.hpp
    include/OptimizerParams/BFGS_Params.hpp

    # Optimizer Base
    include/OptimizerBase/OptimizerBase.hpp
    include/OptimizerBase/LineSearchOptimizer.hpp
    include/OptimizerBase/TrustRegionOptimizer.hpp

    # Line Search Methods
    include/LineSearchMethods/CyclicCoordinateOptimizer.hpp
    include/LineSearchMethods/SteepestDescentOptimizer.hpp
    include/LineSearchMethods/NewtonOptimizer.hpp
    include/LineSearchMethods/ConjuateGradientOptimizer.hpp
    include/LineSearchMethods/MomentumOptimizer.hpp
    include/LineSearchMethods/NesterovMomentumOptimizer.hpp
    include/LineSearchMethods/AdagradOptimizer.hpp
    include/LineSearchMethods/RMSPropOptimizer.hpp
    include/LineSearchMethods/AdaDeltaOptimizer.hpp
    include/LineSearchMethods/AdamOptimizer.hpp
    include/LineSearchMethods/HookeJeevesOptimizer.hpp

    # Newton's method and Quasi-Newton method
    include/LineSearchMethods/NewtonOptimizer.hpp
    include/LineSearchMethods/DFP_Optimizer.hpp
    include/LineSearchMethods/BFGS_Optimizer.hpp

    # Trust Region Methods
    include/TrustRegionMethods/LevenbergMarquardtOptimizer.hpp

    # One Dimensional Search Methods
    include/StepsizeSearch/Accurate/AccurateSearchBase.hpp
    include/StepsizeSearch/Inaccurate/InaccurateSearchBase.hpp

    # Stepsize search base
    include/StepsizeSearch/StepsizeSearchBase.hpp

    # Accurate stepsize search
    include/StepsizeSearch/Accurate/GoldenSectionMethod.hpp
    include/StepsizeSearch/Accurate/FibonacciMethod.hpp
    include/StepsizeSearch/Accurate/BisectionMethod.hpp
    include/StepsizeSearch/Accurate/DichotomousMethod.hpp
    include/StepsizeSearch/Accurate/NewtonMethod.hpp

    # Inaccurate stepsize search
    include/StepsizeSearch/Inaccurate/InaccurateSearchBase.hpp
    include/StepsizeSearch/Inaccurate/GoldsteinMethod.hpp
    include/StepsizeSearch/Inaccurate/ArmijoMethod.hpp
    include/StepsizeSearch/Inaccurate/WolfePowellMethod.hpp
    )

# Examples
#if(CXX_BUILD_EXAMPLES)
#  add_executable(example_robot1 "examples/Robot1/main.cpp")
#endif()

add_executable(${PROJECT_NAME}
    "src/main.cpp")
