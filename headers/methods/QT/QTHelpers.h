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
#include "Eigen/Eigenvalues"
#include "Eigen/SparseCore"

#include "defines.h"
#include "helpers.h"

namespace QT::hlp {

    ///// typedef /////
    typedef Eigen::SparseMatrix<std::complex<double>> matrixType;

    ///// random vectors /////

    // generates normalized random Eigen::Vectors of a given size
    Eigen::VectorXcd getVector(int size);
/*
    // generates normalized random Eigen::Vectors for each matrix block
    std::vector<Eigen::VectorXcd> getVector(const int &N, const int &SIZE, const std::vector<matrixDataMomentumType> &matrixBlocks);
*/
    // generates normalized random Eigen::Vectors for each matrix block
    std::vector<Eigen::VectorXcd> getVector(const std::vector<matrixType> &matrixBlocks);

    /////////////////////////////// Runge-Kutta 4 ///////////////////////////////

    // does one step of the RunkeKutt-Iteration on a given matrix block with a given vector
    Eigen::VectorXcd rungeKutta4Block(const Eigen::VectorXcd &vec, const Eigen::MatrixXcd &H, const double &step);

    std::vector<Eigen::VectorXcd> normalizedVectorList(const std::vector<Eigen::VectorXcd> &vectors, const int &SIZE);

    ///// C /////

    // fourth order Runge-Kutta to calculate beta dependency of the specific heat
    std::vector<double> rungeKutta4_C(const double &start, const double &end, const double &step, const int &N, const std::vector<matrixType> &matrixList);

    ///// X /////

    // fourth order Runge-Kutta to calculate beta dependency of the susceptibility
    std::vector<double> rungeKutta4_X(const double &start, const double &end, const double &step, const int &N,
                                      const std::vector<matrixType> &H_List, const std::vector<matrixType> &S2_List);

    ///// C and X /////

    // fourth order Runge-Kutta to calculate beta dependency of the specific heat and the susceptibility
    std::vector<std::tuple<double, double>> rungeKutta4_CX(const double &start, const double &end, const double &step, const int &N,
                                                           const std::vector<matrixType> &H_List, const std::vector<matrixType> &S2_List);

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
