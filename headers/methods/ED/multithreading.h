#include <iostream>
#include <vector>
#include <complex>
#include "Eigen/Eigenvalues"
#include <omp.h>
#include <mutex>
#include <thread>
#include <chrono>

#include "methods/EDMethods.h"
#include "EDHelpers.h"
#include "helpers.h"
#include "defines.h"

/////////////////////////////// multi-threading ///////////////////////////////

namespace ED::multi {

    static int CURRENT = 1;
    const int PROGRESSBAR_SEGMENTS = 50;

    static std::mutex coutMutex;
    static std::mutex nextJMutex;

    /////////////////////////////// Delta E, C(J) ///////////////////////////////

    void get_DeltaE_CT_const_momentum(double J, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                      double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                      const double &h, const int &N, const int &SIZE);
/*
    void get_DeltaE_CT_const_parity(double J, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                    double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                    const int &N, const int &SIZE);
*/
    void get_DeltaE_CT_const_SI(double J, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                const double &h, const int &N, const int &SIZE);



    void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END, const double &h,
                               int cores, const double &T, const int &N, const int &SIZE);

    /////////////////////////////// Chi(J) ///////////////////////////////

    void get_XT_const_magnetization(double J, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                                    const int &N, const int &SIZE, const double &T);

    void get_XT_const_momentum(double J, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                               const int &N, const int &SIZE, const double &T);

    void start_XT_const(const int &COUNT, const double &START, const double &END,
                        int &cores, const double &T, const int &N, const int &SIZE);

    /////////////////////////////// spin gap ///////////////////////////////

    void get_SpinGap_momentum(std::vector<std::tuple<double, double>> *outDataSpinGap,
                              const double &J, const int &N, const int &SIZE);

    void get_SpinGap_momentum_with_index(std::vector<std::tuple<double, double, int, int, int, int>> *outDataSpinGap,
                                         const double &J, const int &N, const int &SIZE);

    void get_SpinGap_SI(std::vector<std::tuple<double, double>> *outDataSpinGap,
                        const double &J, const int &N, const int &SIZE);

    void start_SpinGap(const int &COUNT, const double &START, const double &END,
                       int &cores, const int &N, const int &SIZE);

    void start_SpinGap_with_index(const int &COUNT, const double &START, const double &END,
                                  int &cores, const int &N, const int &SIZE);


    void startSusceptibilityMultiJ(const double &J_START, const double &J_END, const int &J_COUNT,
                                   const double &BETA_START, const double &BETA_END, const double &BETA_COUNT,
                                   const int &N, const int &SIZE);

    /////////////////////////////// excitation energies ///////////////////////////////

    void startSpecificHeatMultiJ(const double &J_START, const double &J_END, const int &J_COUNT,
                                 const double &BETA_START, const double &BETA_END, const double &BETA_COUNT,
                                 const int &N, const int &SIZE, const double &h);

}