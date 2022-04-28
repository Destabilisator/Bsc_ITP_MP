cmake_minimum_required(VERSION 3.21)
project(Bsc_ITP_MX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc -O3")
#set(CMAKE_EXE_LINKER_FLAGS -static)

set(EIGEN_DIR "./Eigen/eigen-3.4.0")
include_directories(${EIGEN_DIR})
set(HEADER_DIR "./headers")
include_directories(${HEADER_DIR})

include_directories("./src")

add_executable(Bsc_ITP_MX main.h main.cpp
        headers/helpers.h
        headers/helpers.cpp)
