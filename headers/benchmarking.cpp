#include "benchmarking.h"

#define QT_BENCH
#define ED_BENCH_FITS
#define ED_BENCH_EV
#define ED_BENCH_EV_SG
#define ONLY_MAG_ZERO_BLOCK

namespace bench {

    /////////////////////////////// benchmarking ///////////////////////////////

    void bench_ED_QT_SG_runtime(int N_start, int N_end) {

        auto start_timer = std::chrono::steady_clock::now();

        double J_START = 0.0, J_END = 2.0;
        int J_COUNT = 1;
        double T_START = 0.0, T_END = 50.0;

        for (int N = N_start; N <= N_end; N += 2) {
            int SIZE = (int) std::pow(2, N);
            for (int cores: {1, 2, 5, 10}) { /// 1, 2, 5, 10
                if (N > 20 && cores != 1) { continue;}
                for (double stepsize: {0.1, 0.01}) {
                    if (N > 20 && stepsize < 0.05) { continue;}
                    int T_COUNT = (int) ((T_END - T_START) / stepsize);
                    for (int SAMPLES: {1, 2, 3}) { /// 1, 2, 3
                        if (N > 20 && SAMPLES != 1) { continue;}

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
                    /*
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
                            */
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

    void bench_ED_QT_SG_runtime_mag_zero_block(int N_start, int N_end) {

        auto start_timer = std::chrono::steady_clock::now();

        double J_START = 0.0, J_END = 2.0;
        int J_COUNT = 1;
        double T_START = 0.0, T_END = 50.0;

        for (int N = N_start; N <= N_end; N += 2) {
            int SIZE = (int) std::pow(2, N);
            for (int cores: {1, 2, 5, 10}) { /// 1, 2, 5, 10
                if (N > 20 && cores != 1) { continue;}
                for (double stepsize: {0.1, 0.01}) {
                    if (N > 20 && stepsize < 0.05) { continue;}
                    int T_COUNT = (int) ((T_END - T_START) / stepsize);
                    for (int SAMPLES: {1, 2, 3}) { /// 1, 2, 3
                        if (N > 20 && SAMPLES != 1) { continue;}

                        // QT
#ifdef QT_BENCH
                        std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << "), stepsize = " << stepsize
                                  << ", SAMPLES = " << SAMPLES << ", cores = " << cores;
                        auto start_timer_QT = std::chrono::steady_clock::now();
#ifndef ONLY_MAG_ZERO_BLOCK
                        QT::MS::start_calc_spin_gap(J_START, J_END, J_COUNT, T_START, T_END, stepsize, N, SIZE, SAMPLES, cores);
#else
                        start_calc_spin_gap_mag_zero_block(J_START, J_END, J_COUNT, T_START, T_END, stepsize, N, SIZE, SAMPLES, cores);
#endif

                        auto end_timer_QT = std::chrono::steady_clock::now();
                        std::chrono::duration<double> elapsed_seconds_QT = end_timer_QT - start_timer_QT;
#ifndef ONLY_MAG_ZERO_BLOCK
                        save_bench_val(
                                "runtime/data/QT_SG_step_" + std::to_string(stepsize) + "_SAMPLES_" + std::to_string(SAMPLES) +
                                "_cores_" + std::to_string(cores) + ".txt",
                                std::to_string(N) + "\t" + std::to_string(elapsed_seconds_QT.count()));
#else
                        save_bench_val(
                                "runtime/data/QT_SG_zero_block_step_" + std::to_string(stepsize) + "_SAMPLES_" + std::to_string(SAMPLES) +
                                "_cores_" + std::to_string(cores) + ".txt",
                                std::to_string(N) + "\t" + std::to_string(elapsed_seconds_QT.count()));
#endif
#endif
                    }
                    // ED with fit
#ifdef ED_BENCH_FITS
#ifndef ONLY_MAG_ZERO_BLOCK
                    std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << "), stepsize = " << stepsize << ", cores = " << cores;
                    auto start_timer_ED_MJ = std::chrono::steady_clock::now();

                    ED::multi::startSusceptibilityMultiJ(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, N, SIZE,
                                                         cores);

                    auto end_timer_ED_MJ = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed_seconds_ED_MJ = end_timer_ED_MJ - start_timer_ED_MJ;
                    save_bench_val(
                            "runtime/data/ED_MJ_step_" + std::to_string(stepsize) + "_cores_" + std::to_string(cores) + ".txt",
                            std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_MJ.count()));
#endif
#endif
                }
                // ED from EV
#ifdef ED_BENCH_EV_SG
#ifndef ONLY_MAG_ZERO_BLOCK
                std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << ")" << ", cores = " << cores;
                auto start_timer_ED_SG = std::chrono::steady_clock::now();

                ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);

                auto end_timer_ED_SG = std::chrono::steady_clock::now();
                std::chrono::duration<double> elapsed_seconds_ED_SG = end_timer_ED_SG - start_timer_ED_SG;
                save_bench_val("runtime/data/ED_SG_cores_" + std::to_string(cores) + ".txt",
                               std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_SG.count()));
#endif
#endif
#ifdef ED_BENCH_EV
#ifdef ONLY_MAG_ZERO_BLOCK
                if (cores != 1) { continue;}
                std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << ")" << ", cores = " << cores << " MS ED m_z = 0\n";
                auto start_timer_ED_MS = std::chrono::steady_clock::now();
                ed_ms_getEiValsZeroBlock(1.0, 1.0, N, SIZE);
                auto end_timer_ED_MS = std::chrono::steady_clock::now();
                std::chrono::duration<double> elapsed_seconds_ED_MS = end_timer_ED_MS - start_timer_ED_MS;
                save_bench_val("runtime/data/ED_MS_" + std::to_string(cores) + ".txt",
                               std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_MS.count()));
                if (N%4 == 0) {
                    std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << ")" << ", cores = " << cores << " MS ED m_z = 0\n";
                    auto start_timer_ED_SI = std::chrono::steady_clock::now();
                    ed_si_getEiValsZeroBlock(1.0, 1.0, N, SIZE);
                    auto end_timer_ED_SI = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed_seconds_ED_SI = end_timer_ED_SI - start_timer_ED_SI;
                    save_bench_val("runtime/data/ED_SI_" + std::to_string(cores) + ".txt",
                                   std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_SI.count()));
                }
#endif
#endif
            }
        }

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer - start_timer;
        std::cout << "BENCHMARKING: done\n";
        std::cout << "this took " + formatTime(elapsed_seconds) << "\n";

    }

    void bench_ED_QT_H_S2_memory_usage(int N_start, int N_end) {

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

    /////////////////////////////// saving ///////////////////////////////

    void save_bench_val(const std::string &filename, const std::string &content) {
        std::ofstream outfile;
        outfile.open("./results/benchmarking/" + filename, std::ios_base::app);
        outfile << content << "\n";
        outfile.close();
    }

    /////////////////////////////// functions to benchmark ///////////////////////////////

    void start_calc_spin_gap_mag_zero_block(const double &J_START, const double &J_END, const int &J_COUNT,
                                            const double &BETA_START, const double &BETA_END, const double &BETA_STEP,
                                            const int &N, const int &SIZE, const int &SAMPLES, const int &cores) {

        auto start_timer = std::chrono::steady_clock::now();
        std::cout << "\n" << "spin gap, QT, momentum states, N: " << N << ", step size: " << BETA_STEP << " ..." << std::endl;

        ///// get S2 /////
        std::vector<matrixTypeComplex> S2_List;
        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;
        std::vector<std::vector<std::vector<int>>> states(N + 1, std::vector<std::vector<int>>(N/2));
        std::vector<std::vector<std::vector<int>>> R_vals(N + 1, std::vector<std::vector<int>>(N/2));
        // get states
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
        // loop through indices
        for (int m = 0; m <= N; m++) {
            if (m != N/2) { continue;}
            for (int k = k_lower; k <= k_upper; k++) {
                if (states.at(m).at(k - k_lower).empty()) {continue;}
                std::vector<Trp> S2_MtrxLst = QT::hlp::spinMatrixMomentum(N, k, states.at(m).at(k - k_lower),
                                                                      R_vals.at(m).at(k - k_lower));

                int matSize = (int) R_vals.at(m).at(k - k_lower).size();

                matrixTypeComplex S2_mat(matSize,matSize);
                S2_mat.setFromTriplets(S2_MtrxLst.begin(), S2_MtrxLst.end());
                S2_List.emplace_back(S2_mat);
            }
        }

        ///// gather x-data /////
        std::vector<double> beta_Data;
        double beta = BETA_START - BETA_STEP;
        while (beta <= BETA_END) {
            beta += BETA_STEP;
            beta_Data.emplace_back(beta);
        } beta_Data.shrink_to_fit();

        ///// init J vals /////
#ifdef SG_EE_EVEN_J_DIST
        std::vector<double> J_vals;
        for (int J_pos = 0; J_pos < J_COUNT; J_pos++) {
            J_vals.emplace_back(J_START + (J_END - J_START) * J_pos / J_COUNT);
        } J_vals.shrink_to_fit();
#else
        std::vector<double> J_vals;
        double J_init = J_START;
        while (J_init <= J_END) {
            J_vals.emplace_back(J_init);
//            std::cout << "pushing back J = " << J_init << "\n";
            if (J_init >= 0.6 && J_init < 1.25) {
                J_init += 0.04;
            } else {J_init += 0.1;}
        } J_vals.shrink_to_fit();
#endif

        // init progressbar
        int prgbar_segm = 50;
        int curr = 0;
        coutMutex.lock();
        int _p = (int) ( (float) curr / (float) J_vals.size() * (float) prgbar_segm);
        std::cout << "\r[";
        for (int _ = 0; _ < _p; _++) {
            std::cout << "#";
        } for (int _ = _p; _ < prgbar_segm; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) curr / (float) J_vals.size() * 100.0 ) << "% (" << curr << "/" << J_vals.size() << "), J1/J2 = " << J_vals.at(0) << "     ";
        std::cout.flush();
        curr++;
        coutMutex.unlock();



#pragma omp parallel for default(none) num_threads(cores) shared(J_COUNT, J_START, J_END, N, SIZE, SAMPLES, coutMutex, BETA_START, BETA_END, BETA_STEP, S2_List, beta_Data, curr, prgbar_segm, std::cout, J_vals)
        for (int J_pos = 0; J_pos < J_vals.size(); J_pos++) {
            double J = J_vals.at(J_pos); //J_START + (J_END - J_START) * J_pos / J_COUNT;
            std::vector<matrixTypeComplex> H_List = QT::MS::getHamilton_mag_zero_block(J, 1.0, 0.0, N, SIZE);

            for (int s = 1; s <= SAMPLES; s++) {
                std::vector<double> rawData = QT::hlp::rungeKutta4_X(BETA_START, BETA_END, BETA_STEP, N, H_List, S2_List);
            }

            // progressbar
            coutMutex.lock();
            int p = (int) ( (float) curr / (float) J_vals.size() * (float) prgbar_segm);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < prgbar_segm; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) curr / (float) J_vals.size() * 100.0 ) << "% (" << curr << "/" << J_vals.size() << "), J1/J2 = " << J_vals.at(J_pos) << "     ";
            std::cout.flush();
            curr++;
            coutMutex.unlock();
        }

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

    }

    void ed_ms_getEiValsZeroBlock(const double &J1, const double &J2, const int &N, const int &SIZE) {

        std::vector<std::vector<std::complex<double>>> eiVals;
        std::vector<Eigen::MatrixXcd> matrixBlockU;
        std::vector<Eigen::MatrixXcd> matrixBlockS2;

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<int> states_m;
        std::vector<std::vector<int>> states(N/2);
        std::vector<std::vector<int>> R_vals(N/2);
        ED::fillStates(&states_m, N/2, N, SIZE);

        for (int s : states_m) {
            for (int k = k_lower; k <= k_upper; k++) {
                int R = ED::checkState(s, k, N);
                if (R >= 0) {
                    states.at(k - k_lower).push_back(s);
                    R_vals.at(k - k_lower).push_back(R);
                }
            }
        }

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                ED::momentumStates::momentumBlockSolver_withMatrix(J1, J2, k, states.at(k - k_lower),
                                                                   R_vals.at(k - k_lower), eiVals,
                                                                   matrixBlockU, matrixBlockS2, N, SIZE);
            }
        }

    }

    void ed_si_getEiValsZeroBlock(const double &J1, const double &J2, const int &N, const int &SIZE) {

        const double h = 0.0;
        const int mag = N/2;
        auto *eiVals = new std::vector<double>;

        std::vector<Eigen::MatrixXd> matrixBlocks;

        std::vector<int> states_m;
        ED::fillStates(&states_m, mag, N, SIZE);

        std::vector<int> states;
        std::vector<int> R_vals;
        std::vector<int> m_vals;
        std::vector<int> n_vals;
        std::vector<int> c_vals;

        int numberOfStates = 0;
        const int k_upper = N/4;

        for (int k = 0; k <= k_upper; k++) {

            if (mag == N/2) {
                for (int z : {-1, 1}) {
                    for (int p: {-1, 1}) {
                        for (int s: states_m) {
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
                                    states.push_back(s);
                                    R_vals.push_back(sigma * R);
                                    m_vals.push_back(m);
                                    n_vals.push_back(n);
                                    c_vals.push_back(c);
                                    numberOfStates++;
                                }
                            }
                        }
                        if (!states.empty()) {
                            ED::spinInversion::SIBlockSolver(J1, J2, h, k, p, z, states, R_vals, m_vals, n_vals, c_vals, eiVals, &matrixBlocks, N);
                        }
                        states.clear();
                        R_vals.clear();
                        m_vals.clear();
                        n_vals.clear();
                        c_vals.clear();
                    }
                }
            }
            else {
                for (int p : {-1, 1}) {
                    for (int s : states_m) {
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
                                states.push_back(s);
                                R_vals.push_back(sigma * R);
                                m_vals.push_back(m);
                                numberOfStates++;
                            }
                        }
                    }
                    if (!states.empty()) {
                        ED::parityStates::parityBlockSolver(J1, J2, k, p, states, R_vals, m_vals, eiVals, &matrixBlocks, N);
                    }
                    states.clear();
                    R_vals.clear();
                    m_vals.clear();
                }
            }


        }

        delete eiVals;

    }

}
