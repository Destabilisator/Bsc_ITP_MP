#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Eigenvalues>
#include <mutex>
#include <thread>

static int CURRENT = 1;
static const int PROGRESSBAR_SEGMENTS = 50;

static std::mutex coutMutex;
static std::mutex nextJMutex;

#define PI  3.14159265358979323846

///// output, turn off during multithreading /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

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

/////////////////////////////// naiver Ansatz ///////////////////////////////
namespace naiv {
    void fillHamilton(double **hamilton, const double &J1,const double &J2, const int &N, const int &SIZE);

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
               const double &END, const int &COUNT, const int &cores);
}

/////////////////////////////// fixed magnetization blocks ///////////////////////////////

namespace magnetizationBlocks {
    void fillHamiltonBlock(const double &J1, const double &J2, const std::vector<int> &states, double **hamiltonBlock,
                           const int &N, const int &SIZE);

    void getEiValsFromBlock(const double &J1, const double &J2, const int &m,
                            std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks,
                            const int &N, const int &SIZE);

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE);

    void getEiValsZeroBlock(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                            Eigen::MatrixXcd &matrixBlockU, std::vector<int> *states, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);

    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT, const int &cores);
}

/////////////////////////////// momentum states ///////////////////////////////

namespace momentumStates {
    void fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                           const std::vector<int> &R_vals, std::complex<double> **hamiltonBlock, const int &N,
                           const int &SIZE);

    void momentumBlockSolver(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                             const std::vector<int> &R_vals, std::vector<std::complex<double>> *HEiValList,
                             std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE);

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);

    void startSpecificHeat(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                           const double &END, const int &COUNT, const int &cores);
}

/////////////////////////////// parity states (unfinished) ///////////////////////////////

namespace parityStates {
    void getEiVals(double J1, double J2, std::vector<std::complex<double>> *eiVals,
                           std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE);

    void start(const double &J1, const double &J2, const int &N, const int &SIZE);
}
