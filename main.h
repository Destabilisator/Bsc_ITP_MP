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

//// ED /////
#include "helpers.h"
#include "methods/EDMethods.h"
#include "methods/ED/multithreading.h"
#include "methods/ED/3DPlots.h"
#include "methods/QTMethods.h"

#define multiCalc

#ifndef multiCalc
    //#define naiverAnsatz
    //#define momentumStateAnsatz
    //#define magnetizationBlocksAnsatz
    //#define parityStateAnsatz
#endif

///// global variables /////
int N = 8;  // has to be at least 6 and even to preserve the periodic boundary conditions of the delta chain
            // has to be divisible by 4 to use parity and spin inversion
int SIZE;
double J1 = 1.0, J2 = 1.0;
double J_START = 0.0, J_END = 2.0;
int J_COUNT = 50;

double step_size = 0.01; // 0.01; 0.01 and 0.1 yield similar results for N = 16 and 0.1 yields good results for N = 18
int SAMPLES = 1;

double h = 0.0;
double h_START = 0.0, h_END = 2.0;
int h_COUNT = 50;

const double T = 1.0;
double T_START = J_START;
double T_END = J_END;
int T_COUNT = J_COUNT;

///// multi threading stuff /////
const int cpu_cnt = (int) std::thread::hardware_concurrency();// / 2;
