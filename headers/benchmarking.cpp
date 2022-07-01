#include "benchmarking.h"

namespace bench {

    /////////////////////////////// benchmarking ///////////////////////////////

    void bench_ED_QT_SG(int N_start, int N_end) {

        auto start_timer = std::chrono::steady_clock::now();

        double J_START = 0.0, J_END = 2.0;
        int J_COUNT = 50;
        double T_START = 0.0, T_END = 50.0;

        for (int N = N_start; N <= N_end; N += 2) {
            int SIZE = (int) std::pow(2, N);
            for (int cores: {1}) { /// 1, 2, 5, 10
                for (double stepsize: {0.1, 0.01}) {
                    int T_COUNT = (int) ((T_END - T_START) / stepsize);
                    for (int SAMPLES: {1}) { /// 1, 2, 3

                        // QT
                        std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << "), stepsize = " << stepsize
                                  << ", SAMPLES = " << SAMPLES << ", cores = " << cores;
                        auto start_timer_QT = std::chrono::steady_clock::now();

                        QT::MS::start_calc_spin_gap(J_START, J_END, J_COUNT, T_START, T_END, stepsize, N, SIZE, SAMPLES,
                                                    cores);

                        auto end_timer_QT = std::chrono::steady_clock::now();
                        std::chrono::duration<double> elapsed_seconds_QT = end_timer_QT - start_timer_QT;
                        save_bench_val(
                                "runtime/data/QT_SG_step_" + std::to_string(stepsize) + "_SAMPLES_" + std::to_string(SAMPLES) +
                                "_cores_" + std::to_string(cores) + ".txt",
                                std::to_string(N) + "\t" + std::to_string(elapsed_seconds_QT.count()));
                    }
                    // ED with fit
                    std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << "), stepsize = " << stepsize << ", cores = " << cores;
                    auto start_timer_ED_MJ = std::chrono::steady_clock::now();

                    ED::multi::startSusceptibilityMultiJ(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, N, SIZE,
                                                         cores);

                    auto end_timer_ED_MJ = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed_seconds_ED_MJ = end_timer_ED_MJ - start_timer_ED_MJ;
                    save_bench_val(
                            "runtime/data/ED_MJ_step_" + std::to_string(stepsize) + "_cores_" + std::to_string(cores) + ".txt",
                            std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_MJ.count()));
                }
                // ED from EV
                std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << ")" << ", cores = " << cores;
                auto start_timer_ED_SG = std::chrono::steady_clock::now();

                ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);

                auto end_timer_ED_SG = std::chrono::steady_clock::now();
                std::chrono::duration<double> elapsed_seconds_ED_SG = end_timer_ED_SG - start_timer_ED_SG;
                save_bench_val("runtime/data/ED_SG_cores_" + std::to_string(cores) + ".txt",
                               std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_SG.count()));
            }
        }

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer - start_timer;
        std::cout << "BENCHMARKING: done\n";
        std::cout << "this took " + formatTime(elapsed_seconds) << "\n";

    }

    void bench_ED_QT_memory_usage(int N_start, int N_end) {

        auto start_timer = std::chrono::steady_clock::now();

        for (int N = N_start; N <= N_end; N += 2) {
            int SIZE = (int) std::pow(2, N);

            std::cout << "\n"; // BENCHMARKING: N = " << N << " (" << SIZE << "), sorting states (MS) \n";

            // MS sort states
            int k_lower = -(N + 2) / 4 + 1;
            int k_upper = N / 4;
            std::vector<std::vector<std::vector<int>>> states(N + 1, std::vector<std::vector<int>>(N/2));
            std::vector<std::vector<std::vector<int>>> R_vals(N + 1, std::vector<std::vector<int>>(N/2));
            for (int s = 0; s < SIZE; s++) {
                int m = ED::bitSum(s, N);
                for (int k = k_lower; k <= k_upper; k++) {
                    int R = ED::checkState(s, k, N);
                    if (R >= 0) {
                        states.at(m).at(k - k_lower).push_back(s);
                        R_vals.at(m).at(k - k_lower).push_back(R);
                    }
                }
            }

            /// QT H & S2

            std::cout << "BENCHMARKING: N = " << N << " (" << SIZE << "), H & S2 (QT) \n";

            unsigned long long H_elements_QT = 0;
            unsigned long long S2_elements_QT = 0;

            for (int m = 0; m <= N; m++) {
                for (int k = k_lower; k <= k_upper; k++) {
                    if (states.at(m).at(k - k_lower).empty()) {continue;}
                    std::vector<Trp> H_MtrxLst = QT::MS::fillHamiltonBlock(1.0, 1.0, 0.0, k,
                                                                           states.at(m).at(k - k_lower),
                                                                           R_vals.at(m).at(k - k_lower), N);
                    std::vector<Trp> S2_MtrxLst = QT::hlp::spinMatrixMomentum(N, k, states.at(m).at(k - k_lower),
                                                                          R_vals.at(m).at(k - k_lower));

                    int matSize = (int) R_vals.at(m).at(k - k_lower).size();

                    matrixTypeComplex H_mat(matSize,matSize);
                    H_mat.setFromTriplets(H_MtrxLst.begin(), H_MtrxLst.end());
                    H_elements_QT += H_mat.nonZeros();

                    matrixTypeComplex S2_mat(matSize,matSize);
                    S2_mat.setFromTriplets(S2_MtrxLst.begin(), S2_MtrxLst.end());
                    S2_elements_QT += S2_mat.nonZeros();
                }
            }

            save_bench_val("memoryusage/data/QT_H.txt", std::to_string(N) + "\t" + std::to_string(H_elements_QT));
            save_bench_val("memoryusage/data/QT_S2.txt", std::to_string(N) + "\t" + std::to_string(S2_elements_QT));

            /// ED H & S2 MS

            std::cout << "BENCHMARKING: N = " << N << " (" << SIZE << "), H & S2 (ED MS) \n";

            unsigned long long H_S2_elements_ED_MS = 0;

            for (int m = 0; m <= N; m++) {
                for (int k = k_lower; k <= k_upper; k++) {
                    if (states.at(m).at(k - k_lower).empty()) {continue;}
                    unsigned long size = states.at(m).at(k - k_lower).size();
                    H_S2_elements_ED_MS += (unsigned long long) std::pow(size, 2);
                }
            }

            save_bench_val("memoryusage/data/ED_H_S2_MS.txt", std::to_string(N) + "\t" + std::to_string(H_S2_elements_ED_MS));

            states.clear();
            R_vals.clear();

            /// ED H SI

            // Bedingung für SI
            if (N%4 != 0) { continue;}

            std::cout << "BENCHMARKING: N = " << N << " (" << SIZE << "), H (ED SI) \n";

            unsigned long long H_elements_ED_SI = 0;

            std::vector<std::vector<int>> states_m(N+1);
            for (int s = 0; s < SIZE; s++) {
                states_m.at(ED::bitSum(s, N)).push_back(s);
            }

            for (int mag = 0; mag <= N; mag++) {
                for (int k = 0; k <= k_upper; k++) {
                    if (mag == N/2) {
                        for (int z : {-1, 1}) {
                            for (int p: {-1, 1}) {
                                int numberOfStates = 0;
                                for (int s: states_m.at(mag)) {
                                    for (int sigma: {-1, 1}) {
                                        int R, n, m, mp, mz, mpz;
                                        ED::checkStateSI(s, R, mp, mz, mpz, k, N);
                                        int c = ED::spinInversion::getClass_set_m_n(m, n, mp, mz, mpz);
                                        if ((k == 0 || k == k_upper) && sigma == -1) {continue;}
                                        if (c == 2 || c == 4 || c == 5) {
                                            double Na = ED::spinInversion::getNa(m, n, R, sigma, p, z, k, c, N);
                                            double Na_inv = ED::spinInversion::getNa(m, n, R, -sigma, p, z, k, c, N);
                                            if (std::abs(Na) < EPSILON) { R = -1;}
                                            if (sigma == -1 && std::abs(Na_inv) > EPSILON) { R = -1;}
                                        } else if (c == 3) {
                                            double val = 1.0 + (double) z * std::cos(4 * PI * (double) k * (double) m / (double) N);
                                            if (std::abs(val) < EPSILON) { R = -1;}
                                        }
                                        if (R > 0) {
                                            numberOfStates++;
                                        }
                                    }
                                }
                                H_elements_ED_SI += (unsigned long long) std::pow(numberOfStates, 2);
                            }
                        }
                    }
                    else {
                        for (int p : {-1, 1}) {
                            int numberOfStates = 0;
                            for (int s : states_m.at(mag)) {
                                for (int sigma : {-1, 1}) {
                                    int R, m;
                                    ED::checkState(s, &R, &m, k, N);
                                    if ((k == 0 || k == k_upper) && sigma == -1) {continue;}
                                    if (m != -1) {
                                        double val = (double) sigma * (double) p * std::cos(4 * PI * (double) k * (double) m / (double) N);
                                        if (std::abs(1.0 + val) < 1e-10) {R = -1;}
                                        if (sigma == -1 && abs(1.0 - val) > 1e-10) {R = -1;}
                                    }
                                    if (R > 0) {
                                        numberOfStates++;
                                    }
                                }
                            }
                            H_elements_ED_SI += (unsigned long long) std::pow(numberOfStates, 2);
                        }
                    }
                }
            }

            save_bench_val("memoryusage/data/ED_H_SI.txt", std::to_string(N) + "\t" + std::to_string(H_elements_ED_SI));

        }

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer - start_timer;
        std::cout << "BENCHMARKING: done\n";
        std::cout << "this took " + formatTime(elapsed_seconds) << "\n";

    }

/*
    void bench_ED_QT_SG(int N, int cores, int stepsize, int SAMPLES, int type) {
        double J_START = 0.0, J_END = 2.0;
        int J_COUNT = 50;
        double T_START = 0.0, T_END = 50.0;
        int SIZE = (int) std::pow(2, N);
        int T_COUNT = (int) ((T_END - T_START) / stepsize);

        if (type == 0) {
            // QT
            std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << "), stepsize = " << stepsize
                      << ", SAMPLES = " << SAMPLES << ", cores = " << cores;
            auto start_timer_QT = std::chrono::steady_clock::now();

            QT::MS::start_calc_spin_gap(J_START, J_END, J_COUNT, T_START, T_END, stepsize, N, SIZE, SAMPLES,
                                        cores);

            auto end_timer_QT = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds_QT = end_timer_QT - start_timer_QT;
            save_bench_val(
                    "QT_SG_step_" + std::to_string(stepsize) + "_SAMPLES_" + std::to_string(SAMPLES) +
                    "_cores_" + std::to_string(cores) + ".txt",
                    std::to_string(N) + "\t" + std::to_string(elapsed_seconds_QT.count()));
        } else if (type == 1) {
            // ED with fit
            std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << "), stepsize = " << stepsize << ", cores = "
                      << cores;
            auto start_timer_ED_MJ = std::chrono::steady_clock::now();

            ED::multi::startSusceptibilityMultiJ(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, N, SIZE,
                                                 cores);

            auto end_timer_ED_MJ = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds_ED_MJ = end_timer_ED_MJ - start_timer_ED_MJ;
            save_bench_val(
                    "ED_MJ_step_" + std::to_string(stepsize) + "_cores_" + std::to_string(cores) + ".txt",
                    std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_MJ.count()));
        } else if (type == 2) {
            // ED from EV
            std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << ")" << ", cores = " << cores;
            auto start_timer_ED_SG = std::chrono::steady_clock::now();

            ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);

            auto end_timer_ED_SG = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds_ED_SG = end_timer_ED_SG - start_timer_ED_SG;
            save_bench_val("ED_SG_cores_" + std::to_string(cores) + ".txt",
                           std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_SG.count()));
        } else {
                std::cout << "no test executed\n";
        }
    }
*/
    void save_bench_val(const std::string &filename, const std::string &content) {
        std::ofstream outfile;
        outfile.open("./results/benchmarking/" + filename, std::ios_base::app);
//        if (outfile.fail()) {
//            outfile.open("./benchmarking/" + filename);
//        }
        outfile << content << "\n";
        outfile.close();
    }

}
