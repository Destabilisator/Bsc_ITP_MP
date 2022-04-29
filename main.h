#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <thread>
#include <mutex>
#include <ctime>
#include <complex>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>
#include "headers/helpers.h"

#define PI  3.14159265358979323846

// multi threading stuff
const auto cpu_cnt = std::thread::hardware_concurrency();
std::mutex EiValWriterMutex;
std::mutex MatrixAppenderMutex;
std::mutex coutMutex;
