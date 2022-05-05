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

#include "headers/helpers.h"
#include "methods.h"

//// methods ////
#define multiCalc

#ifndef multiCalc
    //#define naiverAnsatz
    //#define momentumStateAnsatz
    //#define magnetizationBlocksAnsatz
    //#define parityStateAnsatz
#endif

///// global variables /////
int N = 6; // has to be at least 6 and even to preserve the periodic boundary conditions of the delta chain
int SIZE;
double J1 = 1.0, J2 = 1.0;
double J_START = 0.0;
double J_END = 2.0;
int J_COUNT = 50;
int J_CURRENT = 1;
const int PROGRESSBAR_SEGMENTS = 50;

const double T = 1.0;
double T_START = J_START;
double T_END = J_END;
int T_COUNT = J_COUNT;

///// multi threading stuff /////
const auto cpu_cnt = std::thread::hardware_concurrency();
std::mutex coutMutex;
std::mutex nextJMutex;
