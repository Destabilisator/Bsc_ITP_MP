#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Eigenvalues>

#include "helpers.h"

#define PI  3.14159265358979323846

///// output, turn off during multithreading /////
//#define showMatrix
#define saveMatrix
#define showEigenvalues
#define saveEigenvalues

/////////////////////////////// parity states (unfinished) ///////////////////////////////

namespace parityStates {
    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *eiVals,
                   std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);
}
