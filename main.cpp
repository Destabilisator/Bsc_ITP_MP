#include "main.h"

//#define DEBUG
//#define ED_METHODS
#define CLUSTER

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

//    omp_set_num_threads(cpu_cnt);
    omp_set_num_threads(cores);
    //omp_set_num_threads(OUTERMOST_NESTED_THREADS * OUTER_NESTED_THREADS * INNER_NESTED_THREADS);
    omp_set_nested(1);

    std::cout << "N: " << N << "; size: " << SIZE << "\n";
    /*
    std::cout << "default J1 and J2 for plots with J = const: J1 = " << J1 << " and J2 = " << J2 << " (currently unchangeable)\n";
    std::cout << "default magnetic field for plots with h = const: h = " << h << "\n";
    std::cout << "default step size for RK4: " << step_size << "\n";
    std::cout << "ranges (might not be applicable to all calculations:" << "\n";
    std::cout << "J_START = " << J_START << ", J_END = " << J_END << " and J_COUNT = " << J_COUNT << "\n";
    std::cout << "T_START = " << T_START << ", T_END = " << T_END << " and T_COUNT = " << T_COUNT << "\n";
    std::cout << "h_START = " << h_START << ", h_END = " << h_END << " and h_COUNT = " << h_COUNT << "\n";
      */
    std::cout << "using " << cores << " cores (hardware limit: " << std::thread::hardware_concurrency() << " cores)\n";
    if (noX) {std::cout << "skipping susceptibility plots\n";}

    std::cout.flush();
    //std::cout << std::endl;

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
        // C and X with QT
        QT::MS::start_calculation_C_J_const(T_START, T_END, step_size, J1, J2, h, N, SIZE, OUTER_NESTED_THREADS); /////////////////////// T_END * T_END
        QT::MS::start_calculation_X_J_const(T_START, T_END, step_size, J1, J2, N, SIZE, OUTER_NESTED_THREADS); /////////////////////// T_END * T_END

        // excitation energy \Delta E(J) and specific heat C(J), T = const
        ED::multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END, h, cores, T, N, SIZE); /////////////////////// J_END * J_END

        // specific heat C(T), J = const
        if (N%4 == 0) {
            ED::spinInversion::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END, T_COUNT); /////////////////////// T_END * T_END
        } else {
            ED::momentumStates::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END, T_COUNT); /////////////////////// T_END * T_END
        }

        // spin gap E_gap (J)
        ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE); /////////////////////// J_END / 2.0

        // susceptibility \Chi(J), T = const
        if (!noX) {
            ED::multi::start_XT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
        }

        // susceptibility \Chi(T), J = const
        if (!noX) {
            //ED::magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
            ED::momentumStates::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
        }

//        // dispersion with fixed J
        ED::momentumStates::startDispersionPlot(J1, J2, h, N, SIZE);
    }
#endif
#endif
/////////////////////////////// testing ///////////////////////////////
#ifdef DEBUG

    ED::momentumStates::startSusceptibility(1.0, 1.0, 18, (int) std::pow(2, 18), T_START, T_END, T_COUNT);

//    ED::momentumStates::start(J1, J2, h, N, SIZE);

//    ED::multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
//    ED::multi::start_XT_const(J_COUNT, J_START, J_END, cores, T, N, SIZE);
//    ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);
//    ED::plot3D::start_C(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, cores, N, SIZE);
//    ED::plot3D::start_X(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, cores, N, SIZE);

    double stepsize = step_size; //(T_END - T_START) / (double) T_COUNT; // 0.01 //step_size;

//    QT::MS::start_calculation_C_J_const(T_START, T_END, stepsize, J1, J2, h, N, SIZE, OUTER_NESTED_THREADS);
    T_COUNT =  (int) ( (T_END - T_START) / stepsize );
//    ED::spinInversion::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END, T_COUNT);


    QT::MS::start_calculation_X_J_const(T_START, T_END, stepsize, J1, J2, N, SIZE, OUTER_NESTED_THREADS);
    QT::MB::start_calculation_X_J_const(T_START, T_END, stepsize, J1, J2, N, SIZE, OUTER_NESTED_THREADS);
    ED::magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
    ED::momentumStates::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);

#endif
#endif
#ifdef CLUSTER

///// Abweichungen als Funktion der Systemgröße, Mittlungen und Temperatur (nur Rohdaten) /////

    double stepsize = step_size;//(T_END - T_START) / (double) T_COUNT; // 0.01
    T_START = 0.0; T_END = 50.0;
    T_COUNT =  (int) ( (T_END - T_START) / stepsize );

    SAMPLES = 5;
//    if (N >= 22) {SAMPLES = 1;}
    cores = SAMPLES;
    omp_set_num_threads(cores);

    /// C ///
/*
//    for (double ss : {1.0, 0.5 , 0.1, 0.05, 0.005, 0.001}) {
//        QT::MS::start_calculation_C_J_const(T_START, T_END, ss, J1, J2, h, N, SIZE, SAMPLES);
//    }

//    QT::MS::start_calculation_C_J_const(T_START, T_END, stepsize, J1, J2, h, N, SIZE, SAMPLES);
    T_COUNT =  (int) ( (T_END - T_START) / stepsize );
    if (N%4 == 0) {
        ED::spinInversion::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END, T_COUNT);
    } else {
        ED::momentumStates::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END, T_COUNT);
    }
*/
    /// X ///
/*
//    for (double ss : {1.0, 0.5 , 0.1, 0.05, 0.005, 0.001}) {
//        QT::MS::start_calculation_X_J_const(T_START, T_END, ss, J1, J2, N, SIZE, SAMPLES);
//    }

//    QT::MS::start_calculation_X_J_const(T_START, T_END, stepsize, J1, J2, N, SIZE, SAMPLES);
//    T_COUNT =  (int) ( (T_END - T_START) / stepsize );
    ED::momentumStates::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);
*/
    /// C and X ///
/*
    if (N == 18) {
        // ED
        T_COUNT =  (int) ( (T_END - T_START) / 0.01 );
        if (N%4 == 0) {
            ED::spinInversion::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END, T_COUNT); /////////////////////// T_END * T_END
        } else {
            ED::momentumStates::startSpecificHeat(J1, J2, h, N, SIZE, T_START, T_END, T_COUNT); /////////////////////// T_END * T_END
        }
        ED::momentumStates::startSusceptibility(J1, J2, N, SIZE, T_START, T_END, T_COUNT);

        // combined QT
        for (double ss : {1.0, 0.5 , 0.1, 0.05, 0.005, 0.001}) { // 1.0, 0.5 , 0.1, 0.05, 0.005, 0.001
            QT::MS::start_calculation_CX_J_const(T_START, T_END, ss, J1, J2, N, SIZE, 12);
        }
        QT::MS::start_calculation_CX_J_const(T_START, T_END, stepsize, J1, J2, N, SIZE, 12);
    }
*/
    /// spin gap ///
///*
    QT::MS::start_calc_spin_gap(J_START, J_END, J_COUNT, T_START, T_END, stepsize, N, SIZE, SAMPLES, cores);
    if (N <= 18) {
        ED::multi::startSusceptibilityMultiJ(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, N, SIZE);
        ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);
    }
//*/
    /// excitation energies ///
///*
    //T_END = 2 * T_END;

    QT::MS::start_calc_excitation_energies(J_START, J_END, J_COUNT, T_START, T_END, stepsize, N, SIZE, SAMPLES, cores);
    if (N <= 18) {
        ED::multi::startSpecificHeatMultiJ(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, N, SIZE, h);
        ED::multi::start_DeltaE_CT_const(J_COUNT, J_START, J_END, h, cores, T, N, SIZE);
    }
//*/

#endif

    std::cout << std::endl;

    return 0;
}
