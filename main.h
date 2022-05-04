#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <thread>
#include <mutex>
#include <ctime>
#include <complex>
#include <Eigen/Eigenvalues>
//#include <Eigen/SparseCore>
#include "headers/helpers.h"

#define PI  3.14159265358979323846

///// global variables /////
int N = 6; // has to be at least 6 and even to preserve the periodic boundary conditions of the delta chain
int SIZE;
double J_START = 0;
double J_END = 2;
int J_COUNT = 50;
int J_CURRENT = 1;
int PROGRASSBAR_SEGMENTS = 50;

double BETA_START = J_START;
double BETA_END = J_END;
int BETA_COUNT = J_COUNT;

///// multi threading stuff /////
const auto cpu_cnt = std::thread::hardware_concurrency();
std::mutex EiValWriterMutex;
std::mutex MatrixAppenderMutex;
std::mutex coutMutex;
std::mutex nextJMutex;
