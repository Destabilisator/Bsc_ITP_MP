#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Eigenvalues>

#define PI  3.14159265358979323846

///// output, turn off during multithreading /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

/////////////////////////////// naiver Ansatz ///////////////////////////////
namespace naiv {
    void fillHamilton(double **hamilton, const double &J1,const double &J2, const int &N, const int &SIZE);

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE, const double &BETA_START,
               const double &BETA_END, const int &BETA_COUNT, const int &cores);
}

/////////////////////////////// fixed magnetization blocks ///////////////////////////////

namespace magnetizationBlocks {
    void fillHamiltonBlock(const double &J1, const double &J2, const std::vector<int> &states, double **hamiltonBlock,
                           const int &N, const int &SIZE);

    void getEiValsFromBlock(const double &J1, const double &J2, const int &m,
                            std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks,
                            const int &N, const int &SIZE);

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE);

    void getEiValsZeroBlock(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                            Eigen::MatrixXcd &matrixBlockU, std::vector<int> *states, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);

    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &BETA_START,
                           const double &BETA_END, const int &BETA_COUNT, const int &cores);
}

/////////////////////////////// momentum states ///////////////////////////////

namespace momentumStates {
    void fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                           const std::vector<int> &R_vals, std::complex<double> **hamiltonBlock, const int &N,
                           const int &SIZE);

    void momentumBlockSolver(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                             const std::vector<int> &R_vals, std::vector<std::complex<double>> *HEiValList,
                             std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE);

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);

    void startSpecificHeat(const double &J1, const double &J2, const int &N, const int &SIZE, const double &BETA_START,
                           const double &BETA_END, const int &BETA_COUNT, const int &cores);
}

/////////////////////////////// parity states (unfinished) ///////////////////////////////

namespace parityStates {
    void getEiVals(double J1, double J2, std::vector<std::complex<double>> *eiVals,
                           std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);
}
