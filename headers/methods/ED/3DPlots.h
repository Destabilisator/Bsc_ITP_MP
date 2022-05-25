#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Eigenvalues"
#include <mutex>
#include <thread>
#include <chrono>

#include "methods.h"
#include "helpers.h"
#include "defines.h"

static int CURRENT3D = 1;
static const int PROGRESSBAR_SEGMENTS3D = 50;

//static std::mutex coutMutex;
static std::mutex nextJMutex3D;

/////////////////////////////// 3D Plots ///////////////////////////////

namespace ED::plot3D {

    void get_C(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
                      const int &TCOUNT, const double &TSTART, const double &TEND,
                      const int &N, const int &SIZE);

    void get_C_parity(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
                             const int &TCOUNT, const double &TSTART, const double &TEND,
                             const int &N, const int &SIZE);

    void get_C_SI(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
                         const int &TCOUNT, const double &TSTART, const double &TEND,
                         const int &N, const int &SIZE);

    void start_C(const int &JCOUNT, const double &JSTART, const double &JEND,
                               const int &TCOUNT, const double &TSTART, const double &TEND,
                               int &cores, const int &N, const int &SIZE);

    void get_X(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
               const int &TCOUNT, const double &TSTART, const double &TEND,
               const int &N, const int &SIZE);

    void get_X_momentum(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
               const int &TCOUNT, const double &TSTART, const double &TEND,
               const int &N, const int &SIZE);

    void start_X(const int &JCOUNT, const double &JSTART, const double &JEND,
                 const int &TCOUNT, const double &TSTART, const double &TEND,
                 int &cores, const int &N, const int &SIZE);

}