#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include "Eigen/SparseCore"
#include <mutex>
#include <omp.h>

#include "helpers.h"
#include "methods/ED/EDHelpers.h"
#include "QTHelpers.h"
#include "defines.h"

/////////////////////////////// magnetization blocks ///////////////////////////////

namespace QT::MB {

    static std::mutex coutMutex;

    typedef std::tuple<int, int, Eigen::SparseMatrix<double>> matrixDataMomentumType; // m, matrix
    typedef std::tuple<int, int, std::vector<int>, std::vector<int>> indexStateVectorType; // m, states,
    typedef Eigen::SparseMatrix<double> matrixType;

    std::vector<Eigen::VectorXcd> getVector(const std::vector<matrixType> &matrixBlocks);

    std::vector<matrixType> getHamilton(const double &J1, const double &J2, const int &N, const int &SIZE);

    Eigen::MatrixXd fillHamiltonBlock(const double &J1, const double &J2, const std::vector<int> &states, const int &N);

    std::vector<matrixType> getS2(const double &J1, const double &J2, const int &N, const int &SIZE);

    Eigen::MatrixXd fillS2Block(const int &N, const std::vector<int> &states);

    std::vector<double> rungeKutta4_X(const double &start, const double &end, const double &step, const int &N, const std::vector<matrixType> &H_List, const std::vector<matrixType> &S2_List);

    void start_calculation_X_J_const(const double &start, const double &end, const double &step,
                                     const double &J1, const double &J2, const int &N, const int &SIZE, const int &SAMPLES);

}
