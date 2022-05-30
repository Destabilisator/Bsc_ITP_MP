#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include "Eigen/SparseCore"

#include "helpers.h"
#include "methods/ED/EDHelpers.h"
#include "QTHelpers.h"
#include "defines.h"

/////////////////////////////// momentum states ///////////////////////////////

typedef std::tuple<int, int, Eigen::SparseMatrix<std::complex<double>>> matrixDataMomentumType; // m, k, matrix
typedef std::tuple<int, int, std::vector<int>, std::vector<int>> indexStateVectorType; // m, k, states, R_vals
typedef Eigen::SparseMatrix<std::complex<double>> matrixType;

namespace QT::MS {

    // returns the MS as a complex SparseMatrix with indices (m, k, H_block)
    std::vector<matrixDataMomentumType> getIndexAndHamilton(const double &J1, const double &J2, const int &N, const int &SIZE);

    // returns the MS as a complex SparseMatrix
    std::vector<Eigen::SparseMatrix<std::complex<double>>> getHamilton(const double &J1, const double &J2, const int &N, const int &SIZE);

    // returns a vector with (m, k, states, R_vals)
    std::vector<indexStateVectorType> getIndexAndStates(const int &N, const int &SIZE);

    // fills the M,k blocks of the MS
    Eigen::MatrixXcd fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                                       const std::vector<int> &R_vals, const int &N);

    // generates normalized random Eigen::Vectors of a given size
    Eigen::VectorXcd getVector(int size);

    // generates normalized random Eigen::Vectors for each matrix block
    std::vector<Eigen::VectorXcd> getVector(const int &N, const int &SIZE, const std::vector<matrixDataMomentumType> &matrixBlocks);

    std::vector<Eigen::VectorXcd> getVector(const std::vector<matrixType> &matrixBlocks);
/*
    // fourth order Runge-Kutta to calculate beta dependency
    std::vector<std::tuple<double, double>> rungeKutta4_C(const double &start, const double &end, const double &step,
                                                          const double &J1, const double &J2, const int &N,
                                                          const indexStateVectorType& matrixIndex);
*/
    // fourth order Runge-Kutta to calculate beta dependency
    std::vector<std::tuple<double, double>> rungeKutta4_C(const double &start, const double &end, const double &step,
                                                          const double &J1, const double &J2, const int &N, const int &SIZE,
                                                          const std::vector<matrixType> &matrixList);

    // calculate and save specific heat as a function of temperature (beta)
    void start_calculation_C_J_const(const double &start, const double &end, const double &step,
                                     const double &J1, const double &J2, const int &N, const int &SIZE);

}
