#include "main.h"

#define DEBUG
//#define ED_METHODS

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

    for(int n = 6; n <= 32; n += 2) {
        N = n;
        SIZE = (int) std::pow(2, N);
        std::cout << "N = " << N << ", SIZE = " << SIZE << std::endl;

//        std::cout << "momentumStates" << std::endl;
//
//        auto start1 = std::chrono::steady_clock::now();
//
//        // get full moment matrix
//        std::vector<std::complex<double>> HEiValList;
//        std::vector<Eigen::MatrixXcd> matrixBlocks;
//
//        int k_lower = -(N + 2) / 4 + 1;
//        int k_upper = N / 4;
//
//        std::vector<std::vector<std::vector<int>>> states(N + 1, std::vector<std::vector<int>>(N/2));
//        std::vector<std::vector<std::vector<int>>> R_vals(N + 1, std::vector<std::vector<int>>(N/2));
//
//        for (int s = 0; s < SIZE; s++) {
//            int m = ED::bitSum(s, N);
//            for (int k = k_lower; k <= k_upper; k++) {
//                int R = ED::checkState(s, k, N);
//                if (R >= 0) {
//                    states.at(m).at(k - k_lower).push_back(s);
//                    R_vals.at(m).at(k - k_lower).push_back(R);
//                }
//            }
//        }
//
//        for (int m = 0; m <= N; m++) {
//            for (int k = k_lower; k <= k_upper; k++) {
//                ED::momentumStates::momentumBlockSolver(J1, J2, k, states.at(m).at(k - k_lower), R_vals.at(m).at(k - k_lower), HEiValList,
//                                                        matrixBlocks, N, SIZE);
//            }
//        }

//        int offset_blocks = 0;
//        Eigen::MatrixXcd H_moment_Block = Eigen::MatrixXcd::Zero(SIZE, SIZE);
//        for (const Eigen::MatrixXcd &M: matrixBlocks) {
//            H_moment_Block.block(offset_blocks, offset_blocks, M.rows(), M.cols()) = M;
//            offset_blocks += (int) M.rows();
//        }

//        auto end1 = std::chrono::steady_clock::now();
//        std::chrono::duration<double> elapsed_seconds1 = end1-start1;
//        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds1) << "\n";
//
//
        std::cout << "getting full sparse matrix" << std::endl;
        auto start2 = std::chrono::steady_clock::now();
//        Eigen::SparseMatrix<std::complex<double>> Matrix = QT::MS::getHamilton(J1, J2, N, SIZE);
        std::vector<matrixDataMomentumType> MatrixList = QT::MS::getHamilton(J1, J2, N, SIZE);

        QT::MS::getVector(N, SIZE, MatrixList);

//
//        for (int i = 0; i < matrixBlocks.size(); i++) {
//            if (matrixBlocks.at(i).isApprox(Eigen::MatrixXcd(std::get<2>(MatrixList.at(i))), EPSILON)) {
//                std::cout << "Oui\n";
//            } else {
//                std::cout << "no bueno\n";
//            }
//        }

//
//        if (Matrix.isApprox(H_moment_Block, EPSILON)) {
//            std::cout << "Oui\n";
//        } else {
//            std::cout << "no bueno\n";
//        }

        auto end2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds2 = end2-start2;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds2) << "\n";

        std::cout << std::endl;

    }





#endif

    std::cout << std::endl;

    return 0;
}
