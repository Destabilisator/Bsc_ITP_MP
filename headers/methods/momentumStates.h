#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Eigenvalues>

#include "helpers.h"

#define PI  3.14159265358979323846

///// output, turn off during multithreading /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

/////////////////////////////// momentum states ///////////////////////////////

namespace momentumStates {

    void fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                           const std::vector<int> &R_vals, Eigen::MatrixXcd &hamiltonBlock, const int &N,
                           const int &SIZE);

    void momentumBlockSolver(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                             const std::vector<int> &R_vals, std::vector<std::complex<double>> &HEiValList,
                             std::vector<Eigen::MatrixXcd> &matrixBlocks, const int &N, const int &SIZE);

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> &HEiValList,
                   std::vector<Eigen::MatrixXcd> &matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);

    void startSpecificHeat(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                           const double &END, const int &COUNT);

    void startDispersionPlot(const double &J1, const double &J2, const int &N, const int &SIZE);

    Eigen::MatrixXcd spinMatrix(const int &N, const int &k, const std::vector<int> &states,
                                const std::vector<int> &R_vals);

    void momentumBlockSolver_withMatrix(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                                        const std::vector<int> &R_vals, std::vector<std::vector<std::complex<double>>> &eiVals,
                                        std::vector<Eigen::MatrixXcd> &matrixBlockU, std::vector<Eigen::MatrixXcd> &matrixBlockS2,
                                        const int &N, const int &SIZE);

    void getEiValsZeroBlock(const double &J1, const double &J2, std::vector<std::vector<std::complex<double>>> &eiVals,
                            std::vector<Eigen::MatrixXcd> &matrixBlockU, std::vector<Eigen::MatrixXcd> &matrixBlockS2,
                            const int &N, const int &SIZE);

    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT);

}
