#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Eigenvalues"

#include "EDHelpers.h"
#include "helpers.h"
#include "defines.h"

///// output, turn off during multithreading /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

/////////////////////////////// parity states ///////////////////////////////

namespace ED::parityStates {

    double get_gk(int k, int N);

    double getNa(int sigma, int m, int Ra, int k, int p, int N);

    double helement(int a, int b, int l, int q, int k, int p, const std::vector<int> &R_vals, const std::vector<int> &m_vals, int N);

    void fillHamiltonParityBlock(const double &J1, const double &J2, const int &k, const int &p, const std::vector<int> &states,
                                 const std::vector<int> &R_vals, const std::vector<int> &m_vals,
                                 Eigen::MatrixXd &hamiltonBlock, const int &N);

    void parityBlockSolver(const double &J1, const double &J2, const int &k, const int &p, const std::vector<int> &states,
                           const std::vector<int> &R_vals, const std::vector<int> &m_vals, std::vector<double> *eiVals,
                           std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N);

    void getEiVals(const double &J1, const double &J2, std::vector<double> *eiVals,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);

}
