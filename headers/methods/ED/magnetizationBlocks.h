#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Eigenvalues"

#include "EDHelpers.h"
#include "helpers.h"

///// output, turn off during multithreading /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

/////////////////////////////// fixed magnetization blocks ///////////////////////////////

namespace ED::magnetizationBlocks {
    
    void fillHamiltonBlock(const double &J1, const double &J2, const std::vector<int> &states, double **hamiltonBlock,
                           const int &N);

    void getEiValsFromBlock(const double &J1, const double &J2, const int &m,
                            std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks,
                            const int &N, const int &SIZE);

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE);

    void getEiValsZeroBlock(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                            Eigen::MatrixXcd &matrixBlockU, const std::vector<int> &states, const int &N);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);

    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT);

}
