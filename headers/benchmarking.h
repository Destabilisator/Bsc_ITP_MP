#include <string>
#include <iostream>
#include <fstream>
#include <sys/resource.h>

#include "methods/EDMethods.h"
#include "methods/ED/multithreading.h"
#include "methods/QTMethods.h"

#include "defines.h"

namespace bench {

    typedef Eigen::Triplet<std::complex<double>> Trp;
    typedef Eigen::SparseMatrix<std::complex<double>> matrixTypeComplex;

    /////////////////////////////// benchmarking ///////////////////////////////

    void bench_ED_QT_SG(int N_start, int N_end);

    void bench_ED_QT_memory_usage(int N_start, int N_end);

    void save_bench_val(const std::string &filename, const std::string &content);

}

