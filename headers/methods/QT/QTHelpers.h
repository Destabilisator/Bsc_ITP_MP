#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <complex>
#include <chrono>
#include <iomanip>
#include <random>
#include <omp.h>
#include "Eigen/Eigenvalues"
#include "Eigen/SparseCore"

#include "defines.h"
#include "helpers.h"
#include "methods/ED/EDHelpers.h"

namespace QT::hlp {

    ///// typedef /////
    typedef Eigen::SparseMatrix<std::complex<double>> matrixTypeComplex;
    typedef Eigen::SparseMatrix<double> matrixTypeReal;
    typedef Eigen::Triplet<std::complex<double>> Trp;

    ///// random vectors /////

    // generates normalized random Eigen::Vectors of a given size
    Eigen::VectorXcd getVector(int size);
/*
    // generates normalized random Eigen::Vectors for each matrix block
    std::vector<Eigen::VectorXcd> getVector(const int &N, const int &SIZE, const std::vector<matrixDataMomentumType> &matrixBlocks);
*/
    // generates normalized random Eigen::Vectors for each matrix block
    std::vector<Eigen::VectorXcd> getVector(const std::vector<matrixTypeComplex> &matrixBlocks);

    std::vector<Eigen::VectorXd> getVector(const std::vector<matrixTypeReal> &matrixBlocks);

    /////////////////////////////// Runge-Kutta 4 ///////////////////////////////

    // does one step of the RunkeKutt-Iteration on a given matrix block with a given vector
    Eigen::VectorXcd rungeKutta4Block(const Eigen::VectorXcd &vec, const matrixTypeComplex &H, const double &step);

    Eigen::VectorXd rungeKutta4Block(const Eigen::VectorXd &vec, const matrixTypeReal &H, const double &step);

    std::vector<Eigen::VectorXcd> normalizedVectorList(const std::vector<Eigen::VectorXcd> &vectors, const int &SIZE);

    ///// C /////

    // fourth order Runge-Kutta to calculate beta dependency of the specific heat
    std::vector<double> rungeKutta4_C(const double &start, const double &end, const double &step, const int &N,
                                      const std::vector<matrixTypeComplex> &matrixList);

    ///// X /////

    // fourth order Runge-Kutta to calculate beta dependency of the susceptibility
    std::vector<double> rungeKutta4_X(const double &start, const double &end, const double &step, const int &N,
                                      const std::vector<matrixTypeComplex> &H_List, const std::vector<matrixTypeComplex> &S2_List);

    // returns the Matrix S^2 for a given set of momentum states
    std::vector<Trp> spinMatrixMomentum(const int &N, const int &k, const std::vector<int> &states, const std::vector<int> &R_vals);

    ///// C and X /////

    // fourth order Runge-Kutta to calculate beta dependency of the specific heat and the susceptibility
    std::vector<std::tuple<double, double>> rungeKutta4_CX(const double &start, const double &end, const double &step, const int &N,
                                                           const std::vector<matrixTypeComplex> &H_List, const std::vector<matrixTypeComplex> &S2_List);

    /////////////////////////////// saving data ///////////////////////////////

    void saveOutData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                     const std::vector<std::tuple<double, double>> &outData, const int &N);

    void saveOutData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                     const std::vector<double> &xData, const std::vector<double> &yData, const int &N);

    void saveOutData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                     const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<double> &yErrData, const int &N);

    void saveStatisticsData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                            const std::vector<double> &xData, const std::vector<std::vector<double>> &RAWyData,
                            const int &SAMPLES, const double &step, const int &N);

    void saveAvgData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                     const std::vector<double> &xData, const std::vector<std::vector<double>> &RAWyData, const int &N);

}
