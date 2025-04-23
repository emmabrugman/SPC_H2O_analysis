cmake_minimum_required (VERSION 3.31)

project (spc)

# Compiler requirements
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

# Find Linear Algebra Header Libraries
find_package(Armadillo REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# Find Linear Algebra Runtime Libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Includes, Links and Defines 
include_directories(${ARMADILLO_INCLUDE_DIRS})
link_libraries(${ARMADILLO_LIBRARIES})
link_libraries(Eigen3::Eigen) 
link_libraries(armadillo)

# add subdir 
add_subdirectory(src)

