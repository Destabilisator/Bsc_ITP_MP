#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Eigenvalues"

#include "EDHelpers.h"
#include "helpers.h"
#include "parityStates.h"
#include "defines.h"

///// output, turn off during multithreading /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

/////////////////////////////// spin inversion ///////////////////////////////

namespace ED::spinInversion {
    double get_gk(int k, int N);

    int getClass_set_m_n(int &m, int &n, const int &mp, const int &mz, const int &mpz);

    double getNa(const int &m, const int &n, const unsigned int &R, const int &sigma, const int &p, const int &z, const int &k, const int &c, const int &N);

    double helement(const int &a, const int &b, const int &l, const int &q, const int &g, const int &k, const int &p,
                    const int &z, const std::vector<int> &R_vals, const std::vector<int> &m_vals,
                    const std::vector<int> &n_vals, const std::vector<int> &c_vals, int N);

    void fillHamiltonSIBlock(const double &J1, const double &J2, const double &h, int k, int p, int z, const std::vector<int> &states,
                             const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                             const std::vector<int> &c_vals, Eigen::MatrixXd &hamiltonBlock, const int &N);

    void SIBlockSolver(const double &J1, const double &J2, const double &h, int k, int p, int z, const std::vector<int> &states,
                       const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                       const std::vector<int> &c_vals, std::vector<double> *eiVals, std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N);

    void getEiVals(const double &J1, const double &J2, const double &h, std::vector<double> *eiVals,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE);

    //////////////////////// susceptibility ////////////////////////

    void SIBlockSolver_withMatrix(const double &J1, const double &J2, const double &h, int k, int p, int z, const std::vector<int> &states,
                                  const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                                  const std::vector<int> &c_vals, std::vector<std::vector<double>> &eiVals, std::vector<Eigen::MatrixXd> &matrixUBlocks,
                                  std::vector<Eigen::MatrixXd> &matrixS2Blocks, const int &N);

    void getEiValsZeroBlock(const double &J1, const double &J2, const double &h, std::vector<std::vector<double>> &eiVals, std::vector<Eigen::MatrixXd> &UBlocks,
                            std::vector<Eigen::MatrixXd> &S2Blocks, const int &N, const int &SIZE);

    void getEiValsMagBlock(const double &J1, const double &J2, const double &h, std::vector<double> *eiVals, const int &N, const int &SIZE, const int &mag);

    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT);

    //////////////////////// specific heat ////////////////////////

    void startSpecificHeat(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE, const double &START,
                           const double &END, const int &COUNT);
}
