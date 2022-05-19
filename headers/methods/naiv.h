#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Eigenvalues>

#include "helpers.h"

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

    void start(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
               const double &END, const int &COUNT, const int &cores);

}
