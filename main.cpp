#include "main.h"

//#define DEBUG

int main(int argc, char* argv[]) {

    bool silent = false;
    bool plotsIn3D = false;
    int cores = (int) cpu_cnt;

    validateInput(argc, argv, cpu_cnt, N, SIZE, J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, silent, cores, plotsIn3D, false, J1, J2);

    std::cout << "N: " << N << "\n";
    std::cout << "SIZE: " << SIZE << "\n";
    std::cout << "J_START: " << J_START << "\n";
    std::cout << "J_END: " << J_END << "\n";
    std::cout << "J_COUNT: " << J_COUNT << "\n";
    std::cout << "T_START: " << T_START << "\n";
    std::cout << "T_END: " << T_END << "\n";
    std::cout << "T_COUNT: " << T_COUNT << "\n";
    std::cout << "cores: " << cores << "\n";
    std::cout << "silent: " << silent << "\n";
    std::cout << "plotsIn3D: " << plotsIn3D << "\n";

    return 1;

/////////////////////////////// calculate quantities ///////////////////////////////
#ifndef DEBUG

    if(plotsIn3D) {
        std::cout << "reeeeeeeeeee\n";
        // 3D Plots of specific heat dependent on T and J
        plot3D::start_C(50, 0, 2, 50, 0, 2, cores, N, SIZE);
        //plot3D::start_C(J_COUNT, J_START, J_END, T_COUNT, T_START, T_END, cores, N, SIZE);
        // 3D Plots of susceptibility dependent on T and J
        //plot3D::start_X(J_COUNT, J_START, J_END, T_COUNT, T_START, T_END, cores, N, SIZE);
    } else {
        // excitation energy \Delta E(J) and specific ehat C(J), T = const
        multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
        // susceptibility \Chi(J), T = const
        multi::start_XT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);

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
    }

#endif
/////////////////////////////// testing ///////////////////////////////
#ifdef DEBUG
//    magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    spinInversion::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);

    // individual methods
//    spinInversion::start(J1, J2, N, SIZE);
//    parityStates::start(J1, J2, N, SIZE);
//    momentumStates::start(J1, J2, N, SIZE);
//    magnetizationBlocks::start(J1, J2, N, SIZE);


    // speed test
//    for (int n = 0; n <= 16; n += 1) {
//        int size = (int) std::pow(2, 12);
//        std::cout << "N: " << n << ", size: " << size << "\n";
//        spinInversion::start(J1, J2, 12, size);
//        parityStates::start(J1, J2, 12, size);
//        std::cout << "\n";
//    }
#endif

    std::cout << "\n" << std::endl;

    return 0;
}
