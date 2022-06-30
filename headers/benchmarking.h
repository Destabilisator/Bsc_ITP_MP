#include <string>
#include <iostream>
#include <fstream>

#include "methods/EDMethods.h"
#include "methods/ED/multithreading.h"
#include "methods/QTMethods.h"

#include "defines.h"

namespace bench {

    /////////////////////////////// benchmarking ///////////////////////////////

    void bench_ED_QT_SG(int N_start, int N_end);

    void save_bench_val(const std::string &filename, const std::string &content);

}

