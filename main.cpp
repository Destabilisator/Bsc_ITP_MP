#include "main.h"

#define DEBUG
//#define ED_METHODS

int main(int argc, char* argv[]) {

    bool silent = false;
    bool plotsIn3D = false;
    int cores = (int) cpu_cnt;
    bool noX = false;

    validateInput(argc, argv, cpu_cnt, N, SIZE, J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, silent, cores, plotsIn3D, true, J1, J2, noX);

    omp_set_num_threads((int) cpu_cnt);

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

   double stepsize = (T_END - T_START) / (double) T_COUNT;
    QT::MS::start_calculation_C_J_const(T_START, T_END, stepsize, J1, J2, N, SIZE, 10);
    ED::momentumStates::startSpecificHeat(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    QT::multi::start_C_J_const(T_START, T_END, stepsize, cores, J1 / J2, N, SIZE);

#endif

    std::cout << std::endl;

    return 0;
}
