#规定cmake的最低版本要求
cmake_minimum_required (VERSION 3.10)
#设置c++编译器
set(CMAKE_CXX_COMPILER "g++")
#项目的名称，可以和文件夹名称（HELLO）不同
project(TKDE_MSPA)
# find_package(OpenMP REQUIRED)
#添加头文件的搜索路径
# include_directories(./include)
include_directories(./src)
set(CMAKE_BINARY_DIR ${PROJECT_SOURCE_DIR}/CodeStore/TKDE_MSPA_CODE)
# include_directories(./input)
#将源文件列表写在变量SrcFiles中
# aux_source_directory(./include INCLUDE)
# aux_source_directory(./src SRCS)
# aux_source_directory(./test TESTS)
#设置可执行文件输出路径
# set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/MSPA.out)
# set(CMAKE_BUILD_TYPE "Debug")
add_definitions(-DRECOVER_MSPATREE)
#设置可执行文件的名称，make之后bin目录下出现myprint.exe
add_executable(MSPA ${PROJECT_SOURCE_DIR}/src/MSPA.cpp)
