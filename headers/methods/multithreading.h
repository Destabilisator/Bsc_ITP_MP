#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Eigenvalues"
#include <mutex>
#include <thread>
#include <chrono>

#include "methods.h"
#include "helpers.h"

static int CURRENT = 1;
static const int PROGRESSBAR_SEGMENTS = 50;

static std::mutex coutMutex;
static std::mutex nextJMutex;

#define PI  3.14159265358979323846

/////////////////////////////// multi-threading ///////////////////////////////

namespace multi {
    void get_DeltaE_CT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                             double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                             const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE);

    void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END, const unsigned int *cpu_cnt,
                               int *cores, const double &T, const int &N, const int &SIZE);

    void get_XT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                      const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE, const double &T);

    void start_XT_const(const int &COUNT, const double &START, const double &END, const unsigned int *cpu_cnt,
                        int *cores, const double &T, const int &N, const int &SIZE);

}