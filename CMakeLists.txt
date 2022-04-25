cmake_minimum_required(VERSION 3.21)
project(Bsc_ITP_MX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")
set(CMAKE_EXE_LINKER_FLAGS -static)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(Bsc_ITP_MX main.cpp)

target_link_libraries (example Eigen3::Eigen)