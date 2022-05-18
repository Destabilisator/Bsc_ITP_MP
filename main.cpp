#include "main.h"

//#define DEBUG

int main(int argc, char* argv[]) {

    bool silent = false;
    int cores = (int) cpu_cnt;

    validateInput(argc, argv, N, SIZE, J_START, J_END, J_COUNT, cpu_cnt, silent, cores, J1, J2, true);

    // syncing T and J ranges
    T_START = J_START;
    T_END = J_END;
    T_COUNT = J_COUNT;

/////////////////////////////// calculate quantities ///////////////////////////////
#ifndef DEBUG
    // excitation energy \Delta E(J) and specific ehat C(J), T = const
    multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END, cpu_cnt, cores, T, N, SIZE);
    // susceptibility \Chi(J), T = const
    multi::start_XT_const(J_COUNT, J_START, J_END, cpu_cnt, cores, T, N, SIZE);

    // specific heat C(T), J = const
    if (N % 4 == 0) {
        spinInversion::startSpecificHeat(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
    } else {
        momentumStates::startSpecificHeat(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
    }

    // susceptibility \Chi(T), J = const
    magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    if (N % 4 == 0) {
//        spinInversion::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    } else {
//        magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    }


    // dispersion with fixed J
    momentumStates::startDispersionPlot(J1, J2, N, SIZE);
#endif
/////////////////////////////// testing ///////////////////////////////
#ifdef DEBUG
    magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
    spinInversion::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);

    // individual methods
//    spinInversion::start(J1, J2, N, SIZE);
//    parityStates::start(J1, J2, N, SIZE);
//    momentumStates::start(J1, J2, N, SIZE);
//    magnetizationBlocks::start(J1, J2, N, SIZE);


    // speed test
//    for (int n = 8; n <= 16; n += 4) {
//        int size = (int) std::pow(2, n);
//        std::cout << "N: " << n << ", size: " << size << "\n";
//        spinInversion::start(J1, J2, n, size);
//        parityStates::start(J1, J2, n, size);
//        std::cout << "\n";
//    }
#endif

    std::cout << "\n" << std::endl;

    return 0;
}
