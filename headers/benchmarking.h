#include <string>
#include <iostream>
#include <fstream>
#include <sys/resource.h>

#include "methods/EDMethods.h"
#include "methods/ED/multithreading.h"
#include "methods/QTMethods.h"

#include "defines.h"

namespace bench {

    static std::mutex coutMutex;
    typedef Eigen::Triplet<std::complex<double>> Trp;
    typedef Eigen::SparseMatrix<std::complex<double>> matrixTypeComplex;

    /////////////////////////////// benchmarking ///////////////////////////////

    void bench_ED_QT_SG_runtime(int N_start, int N_end, int runs);

    void bench_ED_QT_SG_runtime_mag_zero_block(int N_start, int N_end, int runs);

    void bench_ED_QT_SG_runtime_mag_zero_block_quick(int N_start, int N_end, int runs);

    void bench_ED_QT_H_S2_memory_usage(int N_start, int N_end);

    /////////////////////////////// saving ///////////////////////////////

    void save_bench_val(const std::string &filename, const std::string &content);

    /////////////////////////////// functions to benchmark ///////////////////////////////

    void start_calc_spin_gap_mag_zero_block(const double &J_START, const double &J_END, const int &J_COUNT,
                                            const double &BETA_START, const double &BETA_END, const double &BETA_STEP,
                                            const int &N, const int &SIZE, const int &SAMPLES, const int &cores);

    void ed_ms_getEiValsZeroBlock(const double &J1, const double &J2, const int &N, const int &SIZE);

    void ed_si_getEiValsZeroBlock(const double &J1, const double &J2, const int &N, const int &SIZE);

}

