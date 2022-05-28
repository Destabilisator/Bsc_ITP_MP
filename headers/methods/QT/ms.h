#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include "Eigen/SparseCore"

#include "methods/ED/EDHelpers.h"
#include "helpers.h"
#include "defines.h"

/////////////////////////////// momentum states ///////////////////////////////

typedef std::tuple<int, int, Eigen::SparseMatrix<std::complex<double>>> matrixDataMomentumType; // m, k, matrix

namespace QT::MS {

    // returns the MS as a complex SparseMatrix
    std::vector<matrixDataMomentumType> getHamilton(const double &J1, const double &J2, const int &N, const int &SIZE);

    // fills the M,k blocks of the MS
    Eigen::MatrixXcd fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                                       const std::vector<int> &R_vals, const int &N, const int &SIZE);

    // generates normalized random Eigen::Vectors for each matrix block
    std::vector<Eigen::VectorXcd> getVector(const int &N, const int &SIZE, const std::vector<matrixDataMomentumType> &matrixBlocks);

}
