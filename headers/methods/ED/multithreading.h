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

static int CURRENT = 1;
static const int PROGRESSBAR_SEGMENTS = 50;

//static std::mutex coutMutex;
static std::mutex nextJMutex;

/////////////////////////////// multi-threading ///////////////////////////////

namespace ED::multi {

    void get_DeltaE_CT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                             double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                             const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE);

    void get_DeltaE_CT_const_parity(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                    double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                    const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE);
    void get_DeltaE_CT_const_SI(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE);

    void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END,
                               int &cores, const double &T, const int &N, const int &SIZE);

    void get_XT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                      const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE, const double &T);

    void get_XT_const_momentum(double J, int pos, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                               const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE, const double &T);

    void start_XT_const(const int &COUNT, const double &START, const double &END,
                        int &cores, const double &T, const int &N, const int &SIZE);

    void get_SpinGap_momentum(double J, int pos, std::vector<std::tuple<double, double>> *outDataSpinGap,
                              const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE);

    void get_SpinGap_momentum_with_index(double J, int pos, std::vector<std::tuple<double, double, int, int, int, int>> *outDataSpinGap,
                                     const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE);

    void get_SpinGap_SI(double J, int pos, std::vector<std::tuple<double, double>> *outDataSpinGap,
                        const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE);

    void start_SpinGap(const int &COUNT, const double &START, const double &END,
                       int &cores, const int &N, const int &SIZE);

    void start_SpinGap_with_index(const int &COUNT, const double &START, const double &END,
                              int &cores, const int &N, const int &SIZE);

}