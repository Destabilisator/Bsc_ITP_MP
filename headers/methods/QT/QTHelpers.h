#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <complex>
#include <chrono>
#include <iomanip>
#include "Eigen/Eigenvalues"

#include "defines.h"

namespace QT::hlp {

    /////////////////////////////// Runge-Kutta 4 ///////////////////////////////

    Eigen::VectorXcd rungeKutta4Block(const Eigen::VectorXcd &vec, const Eigen::MatrixXcd &H, const double &step);

    std::vector<Eigen::VectorXcd> normalizedVectorList(const std::vector<Eigen::VectorXcd> &vectors, const int &SIZE);

    /////////////////////////////// saving data ///////////////////////////////

    void saveOutData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                     const std::vector<std::tuple<double, double>> &outData, const int &N);

}