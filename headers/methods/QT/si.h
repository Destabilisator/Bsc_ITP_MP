#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/SparseCore"
#include <omp.h>
#include <mutex>

#include "helpers.h"
#include "methods/ED/EDHelpers.h"
#include "methods/EDMethods.h"
#include "QTHelpers.h"
#include "defines.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::SI {

    ///// typedef and global (within namespace) /////

    static std::mutex coutMutex;
    typedef std::tuple<int, int, Eigen::SparseMatrix<std::complex<double>>> matrixDataMomentumType; // m, k, matrix
    typedef std::tuple<int, int, std::vector<int>, std::vector<int>> indexStateVectorType; // m, k, states, R_vals
    typedef Eigen::SparseMatrix<double> matrixType;

    ///// hamilton /////

    // returns the MS as a complex SparseMatrix
    std::vector<matrixType> getHamilton(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE);

    // fills the M,k,p blocks of the PS
    Eigen::MatrixXd fillHamiltonPSBlock(const double &J1, const double &J2, const double &h, const int &k, const int &p,
                                        const std::vector<int> &states, const std::vector<int> &R_vals,
                                        const std::vector<int> &m_vals, const int &N);

    // fills the M,k,p,z blocks of the SI
    Eigen::MatrixXd fillHamiltonSIBlock(const double &J1, const double &J2, const double &h, const int &k, const int &p,
                                        const int &z, const std::vector<int> &states, const std::vector<int> &R_vals,
                                        const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                                        const std::vector<int> &c_vals, const int &N);

    ///// C /////
/*
    // calculate and save the specific heat as a function of temperature (beta)
    void start_calculation_C_J_const(const double &start, const double &end, const double &step, const double &J1,
                                     const double &J2, const double &h, const int &N, const int &SIZE, const int &SAMPLES);

    ///// X /////

    // returns S2 in the MS-Basis as a complex SparseMatrix
    std::vector<matrixType> getS2(const double &J1, const double &J2, const int &N, const int &SIZE);

    Eigen::MatrixXcd fillS2Block(const int &k, const std::vector<int> &states, const std::vector<int> &R_vals, const int &N);

    // calculate and save the susceptibility as a function of temperature (beta)
    void start_calculation_X_J_const(const double &start, const double &end, const double &step,
                                     const double &J1, const double &J2, const int &N, const int &SIZE, const int &SAMPLES);

    ///// C and X /////

    // calculate and save the specific heat and the susceptibility as a function of temperature (beta) with h = 0.0
    void start_calculation_CX_J_const(const double &start, const double &end, const double &step, const double &J1,
                                      const double &J2, const int &N, const int &SIZE, const int &SAMPLES);

    ///// spin gap /////

    // calculate and save the susceptibility as a function of temperature (beta) for multiple J1/J2
    void start_calc_spin_gap(const double &J_START, const double &J_END, const int &J_COUNT,
                             const double &BETA_START, const double &BETA_END, const double &BETA_STEP,
                             const int &N, const int &SIZE, const int &SAMPLES);
*/
}
