#include "main.h"

//#define DEBUG
#define ED_METHODS

int main(int argc, char* argv[]) {

    bool silent = false;
    bool plotsIn3D = false;
    int cores = (int) cpu_cnt;
    bool noX = false;

    validateInput(argc, argv, cpu_cnt, N, SIZE, J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, silent, cores, plotsIn3D, true, J1, J2, noX);

/////////////////////////////// calculate quantities ///////////////////////////////
#ifndef DEBUG
#ifdef ED_METHODS
    if(plotsIn3D) {
        // 3D Plots of specific heat dependent on T and J
        ED::plot3D::start_C(J_COUNT, J_START, J_END, T_COUNT, T_START, T_END, cores, N, SIZE);
        // 3D Plots of susceptibility dependent on T and J
        if (!noX) {
            ED::plot3D::start_X(J_COUNT, J_START, J_END, T_COUNT, T_START, T_END, cores, N, SIZE);
        }
    } else {
        // excitation energy \Delta E(J) and specific heat C(J), T = const
        ED::multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
        // susceptibility \Chi(J), T = const
        if (!noX) {
            ED::multi::start_XT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
        }

        // spin gap E_{gap} (J)
        ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);

        // specific heat C(T), J = const
        if (N % 4 == 0) {
            ED::spinInversion::startSpecificHeat(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
        } else {
            ED::momentumStates::startSpecificHeat(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
        }

        // susceptibility \Chi(T), J = const
        if (!noX) {
            ED::magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
        }
//    if (N % 4 == 0) {
//        ED::spinInversion::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    } else {
//        ED::magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    }

        // dispersion with fixed J
        ED::momentumStates::startDispersionPlot(J1, J2, N, SIZE);
    }
#endif
#endif
/////////////////////////////// testing ///////////////////////////////
#ifdef DEBUG
//    cores = 1;
//    ED::multi::start_SpinGap_with_index(50, 0.01, 2.5, cores, N, SIZE);
    ED::multi::start_SpinGap_with_index(J_COUNT, J_START, J_END, cores, N, SIZE);
//    ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);

    //N = 8; SIZE = (int) std::pow(2,N); T_START = 0; T_END = 20; T_COUNT = 50; J_START = 0; J_END = 2; J_COUNT = 50;
//    ED::multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);

    // susceptibilities
//    ED::magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    ED::momentumStates::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    ED::spinInversion::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);

    // individual methods
//    ED::spinInversion::start(J1, J2, N, SIZE);
//    ED::parityStates::start(J1, J2, N, SIZE);
//    ED::momentumStates::start(J1, J2, N, SIZE);
//    ED::magnetizationBlocks::start(J1, J2, N, SIZE);


    // speed test
//    for (int n = 0; n <= 16; n += 1) {
//        int size = (int) std::pow(2, 12);
//        std::cout << "N: " << n << ", size: " << size << "\n";
//        ED::magnetizationBlocks::startSusceptibility(1.0, 1.0, n, size, 0, 5, 50000);
//        ED::momentumStates::startSusceptibility(1.0, 1.0, n, size, 0, 5, 50000);
////        ED::spinInversion::start(J1, J2, 12, size);
////        ED::parityStates::start(J1, J2, 12, size);
//        std::cout << "\n";
//    }
#endif

    std::cout << std::endl;

    return 0;
}
