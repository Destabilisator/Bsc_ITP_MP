#include "3DPlots.h"

/////////////////////////////// 3D Plots ///////////////////////////////

namespace plot3D {

    void get_C(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
                      const int &TCOUNT, const double &TSTART, const double &TEND,
                      const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex3D.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS3D; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << JSTART << " (" << 0 << "/" << JCOUNT << ")     ";
        std::cout.flush();
        nextJMutex3D.unlock();

        while (true) {

            std::vector<std::complex<double>> eiVals;
            std::vector<Eigen::MatrixXcd> matrixBlocks;

            momentumStates::getEiVals(J, 1.0, &eiVals, &matrixBlocks, N, SIZE);

            // sort eigenvalues
//            eiVals.shrink_to_fit();
//            std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
//                return std::real(c1) < std::real(c2);
//            });

            std::vector<std::tuple<double, double>> C_func_T;

            nextJMutex3D.lock();
            for (int i = 0; i <= TCOUNT; i++) {
                double T = TSTART + (TEND - TSTART) * i / TCOUNT;
                C_func_T.emplace_back(T, getSpecificHeat(T, eiVals, N));
            }

            nextJMutex3D.unlock();

            save3DPlotDataC(J, N, C_func_T);

            // progressbar
            nextJMutex3D.lock();
            int prg = std::min({CURRENT3D, JCOUNT});
            int p = (int) ( (float) prg / (float) JCOUNT * (float) PROGRESSBAR_SEGMENTS3D);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS3D; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) JCOUNT * 100.0 ) << "% J1/J2 = "
                        << JSTART + (JEND - JSTART) * prg / JCOUNT << " (" << prg << "/" << JCOUNT << ")     ";
            std::cout.flush();
            pos = CURRENT3D;
            CURRENT3D++;
            nextJMutex3D.unlock();

            if (pos > JCOUNT) {
                break;
            } else {
                J = JSTART + (JEND - JSTART) * pos / JCOUNT;
            }

        }

    }

    void get_C_parity(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
                             const int &TCOUNT, const double &TSTART, const double &TEND,
                             const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex3D.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS3D; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << JSTART << " (" << 0 << "/" << JCOUNT << ")     ";
        std::cout.flush();
        nextJMutex3D.unlock();

        while (true) {

            std::vector<double> eiVals;
            std::vector<Eigen::MatrixXd> matrixBlocks;

            parityStates::getEiVals(J, 1.0, &eiVals, &matrixBlocks, N, SIZE);


//             //sort eigenvalues
//            eiVals.shrink_to_fit();
//            std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
//                return std::real(c1) < std::real(c2);
//            });

            std::vector<std::tuple<double, double>> C_func_T;

            for (int i = 0; i <= TCOUNT; i++) {
                double T = TSTART + (TEND - TSTART) * i / TCOUNT;
                C_func_T.emplace_back(T, getSpecificHeat(T, eiVals, N));
            }

            save3DPlotDataC(J, N, C_func_T);

            // progressbar
            nextJMutex3D.lock();
            int prg = std::min({CURRENT3D, JCOUNT});
            int p = (int) ( (float) prg / (float) JCOUNT * (float) PROGRESSBAR_SEGMENTS3D);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS3D; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) JCOUNT * 100.0 ) << "% J1/J2 = "
                        << JSTART + (JEND - JSTART) * prg / JCOUNT << " (" << prg << "/" << JCOUNT << ")     ";
            std::cout.flush();
            pos = CURRENT3D;
            CURRENT3D++;
            nextJMutex3D.unlock();

            if (pos > JCOUNT) {
                break;
            } else {
                J = JSTART + (JEND - JSTART) * pos / JCOUNT;
            }

        }

    }

    void get_C_SI(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
                         const int &TCOUNT, const double &TSTART, const double &TEND,
                         const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex3D.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS3D; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << JSTART << " (" << 0 << "/" << JCOUNT << ")     ";
        std::cout.flush();
        nextJMutex3D.unlock();

        while (true) {

            std::vector<double> eiVals;
            std::vector<Eigen::MatrixXd> matrixBlocks;

            spinInversion::getEiVals(J, 1.0, &eiVals, &matrixBlocks, N, SIZE);

//            // sort eigenvalues
//            eiVals.shrink_to_fit();
//            std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
//                return std::real(c1) < std::real(c2);
//            });

            std::vector<std::tuple<double, double>> C_func_T;

            for (int i = 0; i <= TCOUNT; i++) {
                double T = TSTART + (TEND - TSTART) * i / TCOUNT;
                C_func_T.emplace_back(T, getSpecificHeat(T, eiVals, N));
            }

            save3DPlotDataC(J, N, C_func_T);

            // progressbar
            nextJMutex3D.lock();
            int prg = std::min({CURRENT3D, JCOUNT});
            int p = (int) ( (float) prg / (float) JCOUNT * (float) PROGRESSBAR_SEGMENTS3D);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS3D; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) JCOUNT * 100.0 ) << "% J1/J2 = "
                        << JSTART + (JEND - JSTART) * prg / JCOUNT << " (" << prg << "/" << JCOUNT << ")     ";
            std::cout.flush();
            pos = CURRENT3D;
            CURRENT3D++;
            nextJMutex3D.unlock();

            if (pos > JCOUNT) {
                break;
            } else {
                J = JSTART + (JEND - JSTART) * pos / JCOUNT;
            }

        }

    }

    void start_C(const int &JCOUNT, const double &JSTART, const double &JEND,
                               const int &TCOUNT, const double &TSTART, const double &TEND,
                               int &cores, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << JCOUNT << " " << JSTART << " " << JEND << " " << TCOUNT << " " << TSTART << " "  << TEND << "\n";

        std::cout << "\n" << "specific heat (3D): calculating:...";

        if (JCOUNT < cores) {
            cores = JCOUNT;
        }

        std::cout << " (" << cores << ") cores";

        std::thread Threads[cores];

        CURRENT3D = 0 + cores;

        if (N%4 == 0) {
            if (N >= 12) {
                std::cout << ", spin inversion\n";
                for (int i = 0; i < cores; i++) {
                    Threads[i] = std::thread(get_C_SI, JSTART + (JEND - JSTART) * i / JCOUNT, i + 1,
                                             JCOUNT, JSTART, JEND, TCOUNT, TSTART, TEND, N, SIZE);
                }
            } else {
                std::cout << ", parity states\n";
                for (int i = 0; i < cores; i++) {
                    Threads[i] = std::thread(get_C_parity, JSTART + (JEND - JSTART) * i / JCOUNT, i + 1,
                                             JCOUNT, JSTART, JEND, TCOUNT, TSTART, TEND, N, SIZE);
                }
            }

        } else {
            std::cout << ", momentum states\n";
            for (int i = 0; i < cores; i++) {
                Threads[i] = std::thread(get_C, JSTART + (JEND - JSTART) * i / JCOUNT, i + 1,
                                         JCOUNT, JSTART, JEND, TCOUNT, TSTART, TEND, N, SIZE);
            }
        }

        for (int i = 0; i < cores; i++) {
            Threads[i].join();
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";
        std::cout << "\n";

    }

    void get_X(double J, int pos, const int &JCOUNT, const double &JSTART, const double &JEND,
               const int &TCOUNT, const double &TSTART, const double &TEND,
               const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex3D.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS3D; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << JSTART << " (" << 0 << "/" << JCOUNT << ")     ";
        std::cout.flush();
        nextJMutex3D.unlock();

        while (true) {

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

            // progressbar
            nextJMutex3D.lock();
            int prg = std::min({CURRENT3D, JCOUNT});
            int p = (int) ( (float) prg / (float) JCOUNT * (float) PROGRESSBAR_SEGMENTS3D);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS3D; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) JCOUNT * 100.0 ) << "% J1/J2 = "
                        << JSTART + (JEND - JSTART) * prg / JCOUNT << " (" << prg << "/" << JCOUNT << ")     ";
            std::cout.flush();
            pos = CURRENT3D;
            CURRENT3D++;
            nextJMutex3D.unlock();

            if (pos > JCOUNT) {
                break;
            } else {
                J = JSTART + (JEND - JSTART) * pos / JCOUNT;
            }

        }

    }

    void start_X(const int &JCOUNT, const double &JSTART, const double &JEND,
                 const int &TCOUNT, const double &TSTART, const double &TEND,
                 int &cores, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "susceptibility (3D): calculating:...";
        if (JCOUNT < cores) {
            cores = JCOUNT;
        }
        std::cout << " (" << cores << ") cores";

        std::thread Threads[cores];

        CURRENT3D = 1 + cores;

        std::cout << ", magnetization blocks\n";
        for (int i = 0; i < cores; i++) {
            Threads[i] = std::thread(get_X, JSTART + (JEND - JSTART) * i / JCOUNT, i + 1,
                                     JCOUNT, JSTART, JEND, TCOUNT, TSTART, TEND, N, SIZE);
        }

        for (int i = 0; i < cores; i++) {
            Threads[i].join();
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";
        std::cout << "\n";

    }

}