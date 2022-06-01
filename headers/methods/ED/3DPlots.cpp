#include "3DPlots.h"

/////////////////////////////// 3D Plots ///////////////////////////////

namespace ED::plot3D {

    /////////////////////////////// C ///////////////////////////////

    void get_C_momentum(double J, const int &TCOUNT, const double &TSTART, const double &TEND, const int &N, const int &SIZE) {

        std::vector<std::complex<double>> eiVals;
        std::vector<Eigen::MatrixXcd> matrixBlocks;

        momentumStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);

        std::vector<std::tuple<double, double>> C_func_T;

        for (int i = 0; i <= TCOUNT; i++) {
            double T = TSTART + (TEND - TSTART) * i / TCOUNT;
            C_func_T.emplace_back(T, getSpecificHeat(1.0 / T, eiVals, N)); // beta = 1.0 / T => 1.0 / beta = T, undo T => beta conversion
        }

        save3DPlotDataC(J, N, C_func_T);

        // write data
        nextJMutex3D.lock();
//        CURRENT3D++;
        nextJMutex3D.unlock();

    }

    void get_C_parity(double J, const int &TCOUNT, const double &TSTART, const double &TEND, const int &N, const int &SIZE) {

        std::vector<double> eiVals;
        std::vector<Eigen::MatrixXd> matrixBlocks;

        parityStates::getEiVals(J, 1.0, &eiVals, &matrixBlocks, N, SIZE);

        std::vector<std::tuple<double, double>> C_func_T;

        for (int i = 0; i <= TCOUNT; i++) {
            double T = TSTART + (TEND - TSTART) * i / TCOUNT;
            C_func_T.emplace_back(T, getSpecificHeat(1.0 / T, eiVals, N)); // beta = 1.0 / T => 1.0 / beta = T, undo T => beta conversion
        }

        save3DPlotDataC(J, N, C_func_T);

        // write data
        nextJMutex3D.lock();
//        CURRENT3D++;
        nextJMutex3D.unlock();

    }

    void get_C_SI(double J, const int &TCOUNT, const double &TSTART, const double &TEND, const int &N, const int &SIZE) {

        std::vector<double> eiVals;
        std::vector<Eigen::MatrixXd> matrixBlocks;

        spinInversion::getEiVals(J, 1.0, &eiVals, &matrixBlocks, N, SIZE);

        std::vector<std::tuple<double, double>> C_func_T;

        for (int i = 0; i <= TCOUNT; i++) {
            double T = TSTART + (TEND - TSTART) * i / TCOUNT;
            C_func_T.emplace_back(T, getSpecificHeat(1.0 / T, eiVals, N)); // beta = 1.0 / T => 1.0 / beta = T, undo T => beta conversion
        }

        save3DPlotDataC(J, N, C_func_T);

        // write data
        nextJMutex3D.lock();
//        CURRENT3D++;
        nextJMutex3D.unlock();

    }

    void start_C(const double &JSTART, const double &JEND, const int &JCOUNT,
                 const double &TSTART, const double &TEND, const int &TCOUNT,
                 int &cores, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "specific heat (3D): calculating:...";

        std::cout << " (" << cores << ") cores";

        if (N%4 == 0) {
            if (N >= 12) {
                std::cout << ", spin inversion\n";
                int curr = 0;
#pragma omp parallel for default(none) shared(JSTART, JEND, JCOUNT, TSTART, TEND, TCOUNT, N, SIZE, coutMutex, std::cout, curr)
                for (int i = 0; i <= JCOUNT; i++) {
                    double J = JSTART + (JEND - JSTART) * i / JCOUNT;
                    coutMutex.lock();
                    int p = (int) ((float) curr / (float) (JCOUNT) * (float) PROGRESSBAR_SEGMENTS);
                    std::cout << "\r[";
                    for (int _ = 0; _ < p; _++) {
                        std::cout << "#";
                    }
                    for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                        std::cout << ".";
                    }
                    std::cout << "] " << int((float) curr / (float) JCOUNT * 100) << "% J1/J2 = " << J << " (" << curr
                              << "/" << JCOUNT << ")     ";
                    std::cout.flush();
                    curr++;
                    coutMutex.unlock();

                    get_C_SI(J, TCOUNT, TSTART, TEND, N, SIZE);
                }
            } else {
                std::cout << ", parity states\n";
                int curr = 0;
#pragma omp parallel for default(none) shared(JSTART, JEND, JCOUNT, TSTART, TEND, TCOUNT, N, SIZE, coutMutex, std::cout, curr)
                for (int i = 0; i <= JCOUNT; i++) {
                    double J = JSTART + (JEND - JSTART) * i / JCOUNT;
                    coutMutex.lock();
                    int p = (int) ((float) curr / (float) (JCOUNT) * (float) PROGRESSBAR_SEGMENTS);
                    std::cout << "\r[";
                    for (int _ = 0; _ < p; _++) {
                        std::cout << "#";
                    }
                    for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                        std::cout << ".";
                    }
                    std::cout << "] " << int((float) curr / (float) JCOUNT * 100) << "% J1/J2 = " << J << " (" << curr
                              << "/" << JCOUNT << ")     ";
                    std::cout.flush();
                    curr++;
                    coutMutex.unlock();

                    get_C_parity(J, TCOUNT, TSTART, TEND, N, SIZE);
                }
            }
        } else {
            std::cout << ", momentum states\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(JSTART, JEND, JCOUNT, TSTART, TEND, TCOUNT, N, SIZE, coutMutex, std::cout, curr)
            for (int i = 0; i <= JCOUNT; i++) {
                double J = JSTART + (JEND - JSTART) * i / JCOUNT;
                coutMutex.lock();
                int p = (int) ((float) curr / (float) (JCOUNT) * (float) PROGRESSBAR_SEGMENTS);
                std::cout << "\r[";
                for (int _ = 0; _ < p; _++) {
                    std::cout << "#";
                }
                for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                    std::cout << ".";
                }
                std::cout << "] " << int((float) curr / (float) JCOUNT * 100) << "% J1/J2 = " << J << " (" << curr
                          << "/" << JCOUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_C_momentum(J, TCOUNT, TSTART, TEND, N, SIZE);
            }
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";
//        std::cout << "\n";

    }

    /////////////////////////////// Chi ///////////////////////////////

    void get_X_magnetization(double J, const int &N, const int &SIZE,
               const double &JSTART, const double &JEND, const int &JCOUNT,
               const double &TSTART, const double &TEND, const int &TCOUNT) {

        std::vector<int> states;
        std::vector<std::complex<double>> eiVals;

        fillStates(&states, N/2, N, SIZE);
        const int statesCount = (int) states.size();

        Eigen::MatrixXcd matrixBlockU;
        magnetizationBlocks::getEiValsZeroBlock(J, 1.0, &eiVals, matrixBlockU, states, N);

        Eigen::MatrixXd S2 = spinMatrix(N, states);
        Eigen::MatrixXcd Matrix_U_inv_S2_U = Eigen::MatrixXcd::Zero(statesCount, statesCount);
        Matrix_U_inv_S2_U = matrixBlockU.adjoint() * S2 * matrixBlockU;

        std::vector<std::tuple<double, double>> X_func_T;

        for (int i = 0; i <= TCOUNT; i++) {
            double T = TSTART + (TEND - TSTART) * i / TCOUNT;
            X_func_T.emplace_back(T, getSusceptibilityDegeneracy(T, Matrix_U_inv_S2_U, eiVals, N));
        }

        save3DPlotDataX(J, N, X_func_T);

        // write data
        nextJMutex3D.lock();
//        CURRENT3D++;
        nextJMutex3D.unlock();

    }

    void get_X_momentum(double J, const int &N, const int &SIZE,
                        const double &JSTART, const double &JEND, const int &JCOUNT,
                        const double &TSTART, const double &TEND, const int &TCOUNT) {

        std::vector<std::vector<std::complex<double>>> eiVals;
        std::vector<Eigen::MatrixXcd> matrixBlockU;
        std::vector<Eigen::MatrixXcd> matrixBlockS2;
        momentumStates::getEiValsZeroBlock(J, 1.0, eiVals, matrixBlockU, matrixBlockS2, N, SIZE);

        std::vector<std::tuple<double, double>> X_func_T;

        std::vector<Eigen::MatrixXcd> Blocks_U_inv_S2_U;
        for(int i = 0; i < matrixBlockU.size(); i++) {
            Eigen::MatrixXcd M = matrixBlockU.at(i).adjoint() * matrixBlockS2.at(i) * matrixBlockU.at(i);
            Blocks_U_inv_S2_U.push_back(M);
        }

        for (int i = 0; i <= TCOUNT; i++) {
            double current = TSTART + (TEND - TSTART) * i / TCOUNT;
            //current_beta = 1 / current_beta;
            X_func_T.emplace_back(current, getSusceptibilityDegeneracy(current, Blocks_U_inv_S2_U, eiVals, N));
        }

        save3DPlotDataX(J, N, X_func_T);

        // write data
        nextJMutex3D.lock();
//        CURRENT3D++;
        nextJMutex3D.unlock();

    }

    void start_X(const double &JSTART, const double &JEND, const int &JCOUNT,
                 const double &TSTART, const double &TEND, const int &TCOUNT,
                 int &cores, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "susceptibility (3D): calculating:...";

        std::cout << " (" << cores << ") cores";

        if (N >= 10) {
            std::cout << ", momentum states\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(JSTART, JEND, JCOUNT, TSTART, TEND, TCOUNT, N, SIZE, coutMutex, std::cout, curr)
            for (int i = 0; i <= JCOUNT; i++) {
                double J = JSTART + (JEND - JSTART) * i / JCOUNT;
                coutMutex.lock();
                int p = (int) ((float) curr / (float) (JCOUNT) * (float) PROGRESSBAR_SEGMENTS);
                std::cout << "\r[";
                for (int _ = 0; _ < p; _++) {
                    std::cout << "#";
                }
                for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                    std::cout << ".";
                }
                std::cout << "] " << int((float) curr / (float) JCOUNT * 100) << "% J1/J2 = " << J << " (" << curr
                          << "/" << JCOUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_X_momentum(J, N, SIZE, JSTART, JEND, JCOUNT, TSTART, TEND, TCOUNT);
            }
        } else {
            std::cout << ", magnetization blocks\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(JSTART, JEND, JCOUNT, TSTART, TEND, TCOUNT, N, SIZE, coutMutex, std::cout, curr)
            for (int i = 0; i <= JCOUNT; i++) {
                double J = JSTART + (JEND - JSTART) * i / JCOUNT;
                coutMutex.lock();
                int p = (int) ((float) curr / (float) (JCOUNT) * (float) PROGRESSBAR_SEGMENTS);
                std::cout << "\r[";
                for (int _ = 0; _ < p; _++) {
                    std::cout << "#";
                }
                for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                    std::cout << ".";
                }
                std::cout << "] " << int((float) curr / (float) JCOUNT * 100) << "% J1/J2 = " << J << " (" << curr
                          << "/" << JCOUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_X_magnetization(J, N, SIZE, JSTART, JEND, JCOUNT, TSTART, TEND, TCOUNT);
            }
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";
//        std::cout << "\n";

    }

}
