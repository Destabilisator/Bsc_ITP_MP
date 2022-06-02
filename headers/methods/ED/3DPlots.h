#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Eigenvalues"
#include <mutex>
#include <thread>
#include <chrono>

#include "methods/EDMethods.h"
#include "EDHelpers.h"
#include "helpers.h"
#include "defines.h"

/////////////////////////////// 3D Plots ///////////////////////////////

namespace ED::plot3D {

    static int CURRENT3D = 1;
    static const int PROGRESSBAR_SEGMENTS = 50;

    static std::mutex coutMutex;
    static std::mutex nextJMutex3D;

    /////////////////////////////// C ///////////////////////////////

    void get_C_momentum(double &J, const double &h, const int &TCOUNT, const double &TSTART, const double &TEND, const int &N, const int &SIZE);
/*
    void get_C_parity(double J, const int &TCOUNT, const double &TSTART, const double &TEND, const int &N, const int &SIZE);
*/
    void get_C_SI(double J, const double &h, const int &TCOUNT, const double &TSTART, const double &TEND, const int &N, const int &SIZE);

    void start_C(const double &JSTART, const double &JEND, const int &JCOUNT,
                 const double &TSTART, const double &TEND, const int &TCOUNT,
                 const double &h, int &cores, const int &N, const int &SIZE);

    /////////////////////////////// Chi ///////////////////////////////

    void get_X_magnetization(double J, const int &N, const int &SIZE,
                             const double &JSTART, const double &JEND, const int &JCOUNT,
                             const double &TSTART, const double &TEND, const int &TCOUNT);

    void get_X_momentum(double J, const int &N, const int &SIZE,
                        const double &JSTART, const double &JEND, const int &JCOUNT,
                        const double &TSTART, const double &TEND, const int &TCOUNT);

    void start_X(const double &JSTART, const double &JEND, const int &JCOUNT,
                 const double &TSTART, const double &TEND, const int &TCOUNT,
                 int &cores, const int &N, const int &SIZE);

}