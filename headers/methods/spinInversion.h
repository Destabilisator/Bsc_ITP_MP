#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Eigenvalues>

#include "helpers.h"

#define PI  3.14159265358979323846
#define epsilon 1e-10

///// output, turn off during multithreading /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

/////////////////////////////// spin inversion (unfinished) ///////////////////////////////

namespace spinInversion {

    void fillHamiltonSIBlock(const double &J1, const double &J2, int k, int p, int z, const std::vector<int> &states,
                             const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                             const std::vector<int> &c_vals, Eigen::MatrixXd &hamiltonBlock, const int &N);

    void SIBlockSolver(const double &J1, const double &J2, int k, int p, int z, const std::vector<int> &states,
                       const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                       const std::vector<int> &c_vals, std::vector<double> *eiVals, std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N);

    void getEiVals(const double &J1, const double &J2, std::vector<double> *eiVals,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);
}
