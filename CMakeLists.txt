cmake_minimum_required(VERSION 3.21)
project(Bsc_ITP_MX)

set(CMAKE_CXX_STANDARD 20)

find_package(Threads)

if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    set(CMAKE_EXE_LINKER_FLAGS -static)
else()
    file(READ "/etc/issue" ETC_ISSUE)
    string(REGEX MATCH "Debian|Ubuntu" DIST ${ETC_ISSUE})

    if(DIST STREQUAL "Ubuntu")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    else()
        find_package(MKL REQUIRED)
        set(MKLROOT /opt/intel/oneapi/mkl/latest)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -m64 -I${MKLROOT}/include -DEIGEN_USE_MKL_ALL") ##############
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -m64 -I${MKLROOT}/include -DEIGEN_USE_MKL_ALL")
    endif()
endif()

set(EIGEN_DIR "./Eigen/eigen-3.4.0")
include_directories(${EIGEN_DIR})
set(HEADER_DIR "./headers")
include_directories(${HEADER_DIR})
set(MAIN_DIR "./")
include_directories(${MAIN_DIR})

include_directories("./src")

add_executable(Bsc_ITP_MX main.h main.cpp
        headers/helpers.h
        headers/helpers.cpp
        headers/methods.h
        headers/methods/multithreading.h
        headers/methods/multithreading.cpp
        headers/methods/naiv.h
        headers/methods/naiv.cpp
        headers/methods/magnetizationBlocks.h
        headers/methods/magnetizationBlocks.cpp
        headers/methods/momentumStates.h
        headers/methods/momentumStates.cpp
        headers/methods/parityStates.h
        headers/methods/parityStates.cpp
        headers/methods/spinInversion.h
        headers/methods/spinInversion.cpp)

if (WIN32)
    target_link_libraries(Bsc_ITP_MX ${CMAKE_THREAD_LIBS_INIT})
elseif (DIST STREQUAL "Ubuntu")
    target_link_libraries(Bsc_ITP_MX ${CMAKE_THREAD_LIBS_INIT})
else()
    target_link_libraries(Bsc_ITP_MX ${CMAKE_THREAD_LIBS_INIT} -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
            ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl)
endif()
