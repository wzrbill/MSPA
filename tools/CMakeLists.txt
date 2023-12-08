cmake_minimum_required (VERSION 3.10)
set(CMAKE_CXX_COMPILER "g++")
project(TKDE_MSPA)
include_directories(./src)
set(CMAKE_BINARY_DIR ${PROJECT_SOURCE_DIR}/CodeStore/TKDE_MSPA_CODE)
#comment RECOVER_MSPATREE if run on graphs larger than demo graph.
add_definitions(-DRECOVER_MSPATREE)
add_executable(MSPA ${PROJECT_SOURCE_DIR}/src/MSPA.cpp)
