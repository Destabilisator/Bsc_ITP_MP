#include "main.h"

#define DEBUG
//#define ED_METHODS
//#define CLUSTER

int main(int argc, char* argv[]) {

    bool silent = false;
    bool plotsIn3D = false;
    int cores = (int) cpu_cnt;
    bool noX = false;

    validateInput(argc, argv, cpu_cnt, N, SIZE,
                  J_START, J_END, J_COUNT,
                  T_START, T_END, T_COUNT,
                  h, h_START, h_END, h_COUNT,
                  silent, cores, plotsIn3D, true, J1, J2, noX);

    omp_set_num_threads(OUTER_NESTED_THREADS * INNER_NESTED_THREADS);
    omp_set_nested(1);

    std::cout << "N: " << N << "; size: " << SIZE << "\n";
    std::cout << "default J1 and J2 for plots with J = const: J1 = " << J1 << " and J2 = " << J2 << " (currently unchangeable)\n";
    std::cout << "default magnetic field for plots with h = const: h = " << h << "\n";
    std::cout << "default step size for RK4: " << step_size << "\n";
    std::cout << "ranges (might not be applicable to all calculations:" << "\n";
    std::cout << "J_START = " << J_START << ", J_END = " << J_END << " and J_COUNT = " << J_COUNT << "\n";
    std::cout << "T_START = " << T_START << ", T_END = " << T_END << " and T_COUNT = " << T_COUNT << "\n";
    std::cout << "h_START = " << h_START << ", h_END = " << h_END << " and h_COUNT = " << h_COUNT << "\n";
    std::cout << "using " << cores << " cores (hardware limit: " << std::thread::hardware_concurrency() << " cores)\n";
    if (noX) {std::cout << "skipping susceptibility plots\n";}
    std::cout << std::endl;

/////////////////////////////// calculate quantities ///////////////////////////////
#ifndef CLUSTER
#ifndef DEBUG
#ifdef ED_METHODS
    if(plotsIn3D) {
        // 3D Plots of specific heat dependent on T and J
        ED::plot3D::start_C(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, h, cores, N, SIZE);
        // 3D Plots of susceptibility dependent on T and J
        if (!noX) {
            ED::plot3D::start_X(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, cores, N, SIZE);
        }
    } else {
        QT::MS::start_calculation_C_J_const(T_START, T_END * T_END, step_size, J1, J2, h, N, SIZE, 10); ///////////////////////

        // excitation energy \Delta E(J) and specific heat C(J), T = const
        ED::multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END * J_END, h, cores, T, N, SIZE); ///////////////////////

        // specific heat C(T), J = const
        if (N%4 == 0) {
            ED::spinInversion::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END * T_END, T_COUNT); ///////////////////////
        } else {
            ED::momentumStates::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END * T_END, T_COUNT); ///////////////////////
        }

        // spin gap E_gap (J)
        ED::multi::start_SpinGap(J_COUNT, J_START, J_END / 2.0, cores, N, SIZE); ///////////////////////

        // susceptibility \Chi(J), T = const
        if (!noX) {
            ED::multi::start_XT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
        }

        // susceptibility \Chi(T), J = const
        if (!noX) {
            ED::magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
        }

//        // dispersion with fixed J
        ED::momentumStates::startDispersionPlot(J1, J2, h, N, SIZE);
    }
#endif
#endif
/////////////////////////////// testing ///////////////////////////////
#ifdef DEBUG

//    ED::momentumStates::start(J1, J2, h, N, SIZE);

//    ED::multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
//    ED::multi::start_XT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
//    ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);
//    ED::plot3D::start_C(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, cores, N, SIZE);
//    ED::plot3D::start_X(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, cores, N, SIZE);

    double stepsize = step_size;//(T_END - T_START) / (double) T_COUNT; // 0.01
    QT::MS::start_calculation_C_J_const(T_START, T_END, stepsize, J1, J2, h, N, SIZE, 24);
    T_COUNT =  (int) ( (T_END - T_START) / stepsize );
    ED::spinInversion::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END, T_COUNT);


//    QT::MS::start_calculation_X_J_const(T_START, T_END, stepsize, J1, J2, N, SIZE, 5);
//    QT::MB::start_calculation_X_J_const(T_START, T_END, stepsize, J1, J2, N, SIZE, 5);
//    ED::magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
//    ED::momentumStates::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);

#endif
#endif
#ifdef CLUSTER
///// Abweichungen als Funktion der Systemgröße, Mittlungen und Temperatur /////
#endif

    std::cout << std::endl;

    return 0;
}
