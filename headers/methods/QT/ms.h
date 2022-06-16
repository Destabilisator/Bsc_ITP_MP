#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/SparseCore"
#include <omp.h>
#include <mutex>

#include "helpers.h"
#include "methods/ED/EDHelpers.h"
#include "QTHelpers.h"
#include "defines.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::MS {

    ///// typedef and global (within namespace) /////

    static std::mutex coutMutex;
    typedef std::tuple<int, int, Eigen::SparseMatrix<std::complex<double>>> matrixDataMomentumType; // m, k, matrix
    typedef std::tuple<int, int, std::vector<int>, std::vector<int>> indexStateVectorType; // m, k, states, R_vals
    typedef Eigen::SparseMatrix<std::complex<double>> matrixTypeComplex;
    typedef Eigen::VectorXcd vectorTypeComplex;
    typedef Eigen::Triplet<std::complex<double>> Trp;

    ///// hamilton /////
/*
    // returns the MS as a complex SparseMatrix with indices (m, k, H_block)
    std::vector<matrixDataMomentumType> getIndexAndHamilton(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE);
*/
    // returns the MS as a complex SparseMatrix
    std::vector<matrixTypeComplex> getHamilton(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE);
/*
    // returns a vector with (m, k, states, R_vals)
    std::vector<indexStateVectorType> getIndexAndStates(const int &N, const int &SIZE);
*/
    // fills the M,k blocks of the MS
    std::vector<Trp> fillHamiltonBlock(const double &J1, const double &J2, const double &h, const int &k, const std::vector<int> &states,
                                       const std::vector<int> &R_vals, const int &N);

    ///// C /////

    // calculate and save the specific heat as a function of temperature (beta)
    void start_calculation_C_J_const(const double &start, const double &end, const double &step, const double &J1,
                                     const double &J2, const double &h, const int &N, const int &SIZE, const int &SAMPLES);

    ///// X /////

    // returns S2 in the MS-Basis as a complex SparseMatrix
    std::vector<matrixTypeComplex> getS2(const double &J1, const double &J2, const int &N, const int &SIZE);

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

    ///// excitation energies /////

    // calculate and save the specific heat as a function of temperature (beta) for multiple J1/J2
    void start_calc_excitation_energies(const double &J_START, const double &J_END, const int &J_COUNT,
                                        const double &BETA_START, const double &BETA_END, const double &BETA_STEP,
                                        const int &N, const int &SIZE, const int &SAMPLES);

}
