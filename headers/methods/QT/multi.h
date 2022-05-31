#include <iostream>
#include <vector>
#include <complex>
#include <mutex>
#include <thread>
#include <chrono>

#include "methods/QTMethods.h"
#include "QTHelpers.h"
#include "helpers.h"
#include "defines.h"

/////////////////////////////// multi-threading ///////////////////////////////

namespace QT::multi {

    static int NEXT = 0;
    static const int PROGRESSBAR_SEGMENTS = 50;
    static std::mutex nextQTMutex;

    typedef std::tuple<int, int, std::vector<int>, std::vector<int>> indexStateVectorType;

    /////////////////////////////// C(T) ///////////////////////////////

    void start_C_J_const(const double &START, const double &END, const double &STEP, int &cores,
                         const double &J, const int &N, const int &SIZE);

    void get_C_J_const(double J, int pos, std::vector<std::vector<std::tuple<double, double>>> *outData,
                       const std::vector<indexStateVectorType> &indexList, const double &START,
                       const double &END, const double &STEP, const int &N, const int &SIZE);

}