#include <iostream>
#include <vector>
#include <complex>
//#include "Eigen/Eigenvalues"
#include "Eigen/SparseCore"

#include "methods/ED/EDHelpers.h"
#include "helpers.h"
#include "defines.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::hamilton {

    // returns the hamilton as a complex SparseMatrix
    Eigen::SparseMatrix<std::complex<double>> getHamilton(const double &J1, const double &J2, const int &N, const int &SIZE);

    // fills the M,k blocks of the hamilton
    void fillHamiltonBlockMomentum(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                                   const std::vector<int> &R_vals, std::vector<Eigen::MatrixXcd> &hamiltonBlocks,
                                   const int &N, const int &SIZE);

}
