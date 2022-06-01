#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <complex>
#include <chrono>
#include <iomanip>
#include <numeric>
#include "Eigen/Eigenvalues"
#include "defines.h"


/////////////////////////////// bits ///////////////////////////////

// print out first N bits of the integer s
void printBits(int s, int N);

/////////////////////////////// statistics ///////////////////////////////

std::tuple<double, double> get_mean_and_se(const std::vector<double> &data);

/////////////////////////////// others ///////////////////////////////

// check validity of user input
// [executable] N J_START J_END J_COUNT CORES SILENT
void validateInput(int &argc, char *argv[], const unsigned int &cpu_cnt, int &N, int &SIZE, double &J_START,
                   double &J_END,
                   int &J_COUNT, double &T_START, double &T_END, int &T_COUNT, bool &silent, int &cores,
                   bool &plotsIn3D,
                   bool skipSilent, const double &J1, const double &J2, bool &noX);

// converts time in s in hh:mm:ss
std::string formatTime(std::chrono::duration<double> elapsed_seconds);
