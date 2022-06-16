#include "ms.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::MS {

    ///// hamilton /////

    std::vector<indexStateVectorType> getIndexAndStates(const int &N, const int &SIZE) {

        std::vector<std::tuple<int, int, std::vector<int>, std::vector<int>>> indexAndStates;

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<std::vector<std::vector<int>>> states(N + 1, std::vector<std::vector<int>>(N/2));
        std::vector<std::vector<std::vector<int>>> R_vals(N + 1, std::vector<std::vector<int>>(N/2));

        int numberOfStates = 0;

        // get states
        for (int s = 0; s < SIZE; s++) {
            int m = ED::bitSum(s, N);
            for (int k = k_lower; k <= k_upper; k++) {
                int R = ED::checkState(s, k, N);
                if (R >= 0) {
                    states.at(m).at(k - k_lower).push_back(s);
                    R_vals.at(m).at(k - k_lower).push_back(R);
                    numberOfStates++;
                }
            }
        }

        // fill return vector
        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                if (states.at(m).at(k - k_lower).empty()) {continue;}
                indexAndStates.emplace_back(m, k, states.at(m).at(k - k_lower), R_vals.at(m).at(k - k_lower));
            }
        }

        return indexAndStates;

    }

    std::vector<matrixDataMomentumType> getIndexAndHamilton(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE) {

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<std::vector<std::vector<int>>> states(N + 1, std::vector<std::vector<int>>(N/2));
        std::vector<std::vector<std::vector<int>>> R_vals(N + 1, std::vector<std::vector<int>>(N/2));

        int numberOfStates = 0;

        // get states
        for (int s = 0; s < SIZE; s++) {
            int m = ED::bitSum(s, N);
            for (int k = k_lower; k <= k_upper; k++) {
                int R = ED::checkState(s, k, N);
                if (R >= 0) {
                    states.at(m).at(k - k_lower).push_back(s);
                    R_vals.at(m).at(k - k_lower).push_back(R);
                    numberOfStates++;
                }
            }
        }

        std::vector<Eigen::MatrixXcd> matrixBlocks;
        std::vector<matrixDataMomentumType> matList;

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                if (states.at(m).at(k - k_lower).empty()) {continue;}
                Eigen::MatrixXcd Mtrx = fillHamiltonBlock(J1, J2, h, k, states.at(m).at(k - k_lower),
                                                          R_vals.at(m).at(k - k_lower), N);
                matList.emplace_back(m, k, Mtrx.sparseView());
            }
        }

        matList.shrink_to_fit();
        return matList;

    }

    std::vector<matrixTypeComplex> getHamilton(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE) {

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

        std::vector<matrixTypeComplex> matList;

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                if (states.at(m).at(k - k_lower).empty()) {continue;}
                Eigen::MatrixXcd Mtrx = fillHamiltonBlock(J1, J2, h, k, states.at(m).at(k - k_lower),
                                                          R_vals.at(m).at(k - k_lower), N);
                matList.emplace_back(Mtrx.sparseView());
            }
        }

        matList.shrink_to_fit();
        return matList;

    }

    Eigen::MatrixXcd fillHamiltonBlock(const double &J1, const double &J2, const double &h, const int &k, const std::vector<int> &states,
                                       const std::vector<int> &R_vals, const int &N) {

        const int statesCount = (int) states.size();
        Eigen::MatrixXcd hamiltonBlock = Eigen::MatrixXcd::Zero(statesCount,statesCount);

        for (int a = 0; a < statesCount; a++) {
            int s = states.at(a);
            int mag = ED::bitSum(s, N) - N/2;
            hamiltonBlock(a,a) += (double) mag * h;
            for (int n = 0; n < N / 2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    hamiltonBlock(a,a) += std::complex<double>(0.25 * J1, 0.0);
                } else {
                    hamiltonBlock(a,a) -= std::complex<double>(0.25 * J1, 0.0);
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int r = 0, l = 0;
                    ED::representative(d, &r, &l, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                        hamiltonBlock(a,b) += (std::complex<double>) 0.5 * J1 * sqrt((double) R_vals.at(a)
                                                                                     / (double) R_vals.at(b)) * std::exp(numC);
                    }
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    hamiltonBlock(a,a) += std::complex<double>(0.25 * J2, 0.0);
                } else {
                    hamiltonBlock(a,a) -= std::complex<double>(0.25 * J2, 0.0);
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    int r = 0, l = 0;
                    ED::representative(d, &r, &l, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                        hamiltonBlock(a,b) += (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a)
                                                                                     / (double) R_vals.at(b)) * std::exp(numC);
                    }
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    hamiltonBlock(a,a) += std::complex<double>(0.25 * J2, 0.0);
                } else {
                    hamiltonBlock(a,a) -= std::complex<double>(0.25 * J2, 0.0);
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    int r = 0, l = 0;
                    ED::representative(d, &r, &l, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                        hamiltonBlock(a,b) += (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a)
                                                                                     / (double) R_vals.at(b)) * std::exp(numC);
                    }
                }
            }
        }

        return hamiltonBlock;

    }

    ///// C /////

    void start_calculation_C_J_const(const double &start, const double &end, const double &step, const double &J1,
                                     const double &J2, const double &h, const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();
        std::cout << "\n" << "C(T), J = const, QT, momentum states, N: " << N << ", step size: " << step << " ..." << std::endl;

        std::vector<matrixTypeComplex> matrixList = getHamilton(J1, J2, h, N, SIZE);

//        typedef std::vector<std::tuple<double, double>> dataVectorType;
        std::vector<std::vector<double>> outData;

        int prgbar_segm = 50;
        int curr = 0;
        coutMutex.lock();
        int _p = (int) ( (float) curr / (float) SAMPLES * (float) prgbar_segm);
        std::cout << "\r[";
        for (int _ = 0; _ < _p; _++) {
            std::cout << "#";
        } for (int _ = _p; _ < prgbar_segm; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) curr / (float) SAMPLES * 100.0 ) << "% (" << curr << "/" << SAMPLES << ")     ";
        std::cout.flush();
        curr++;
        coutMutex.unlock();
#if OUTER_NESTED_THREADS > 1
        #pragma omp parallel for num_threads(OUTER_NESTED_THREADS) default(none) shared(SAMPLES, coutMutex, curr, prgbar_segm, std::cout, start, end, step, N, matrixList, outData)
#endif
        for (int s = 1; s <= SAMPLES; s++) {
//            int p = (int) ( (float) s / (float) SAMPLES * (float) prgbar_segm);
            std::vector<double> rawData = hlp::rungeKutta4_C(start, end, step, N, matrixList);
//            #pragma omp critical

            coutMutex.lock();
            outData.emplace_back(rawData);
            int p = (int) ( (float) curr / (float) SAMPLES * (float) prgbar_segm);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < prgbar_segm; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) curr / (float) SAMPLES * 100.0 ) << "% (" << curr << "/" << SAMPLES << ")     ";
            std::cout.flush();
            curr++;
            coutMutex.unlock();
        }

        outData.shrink_to_fit();

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";
        std::cout << "preparing and saving data\n";

        ///// avg and saving results (for different sample sizes n) /////

        // gather x-data
        std::vector<double> beta_Data;
        double beta = start - step;
        while (beta <= end) {
            beta += step;
            beta_Data.emplace_back(beta);
        } beta_Data.shrink_to_fit();

        hlp::saveStatisticsData("data_specific_heat_J_const_QT",
                                "N: " + std::to_string(N) + "\n"
                                + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                                + "h: " + std::to_string(h)+ "\n"
                                + "samples: " + std::to_string(SAMPLES) + "\n"
                                + "this took: " + formatTime(elapsed_seconds),
                                "beta in kb / J2", "C in J2", beta_Data, outData, SAMPLES, step, N);

    }

    ///// X /////

    void start_calculation_X_J_const(const double &start, const double &end, const double &step, const double &J1,
                                     const double &J2, const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();
        std::cout << "\n" << "X(T), J = const, QT, momentum states, N: " << N << ", step size: " << step << " ..." << std::endl;

        ///// get H and S2 /////
        std::vector<matrixTypeComplex> H_List;
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
            for (int k = k_lower; k <= k_upper; k++) {
                if (states.at(m).at(k - k_lower).empty()) {continue;}
                Eigen::MatrixXcd H_Mtrx = fillHamiltonBlock(J1, J2, 0.0, k, states.at(m).at(k - k_lower),
                                                          R_vals.at(m).at(k - k_lower), N);
                Eigen::MatrixXcd S2_Mtrx = ED::spinMatrixMomentum(N, k, states.at(m).at(k - k_lower), R_vals.at(m).at(k - k_lower));
                H_List.emplace_back(H_Mtrx.sparseView());
                S2_List.emplace_back(S2_Mtrx.sparseView());
            }
        }

        // init progressbar
        int prgbar_segm = 50;
        int curr = 0;
        coutMutex.lock();
        int _p = (int) ( (float) curr / (float) SAMPLES * (float) prgbar_segm);
        std::cout << "\r[";
        for (int _ = 0; _ < _p; _++) {
            std::cout << "#";
        } for (int _ = _p; _ < prgbar_segm; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) curr / (float) SAMPLES * 100.0 ) << "% (" << curr << "/" << SAMPLES << ")     ";
        std::cout.flush();
        curr++;
        coutMutex.unlock();

        std::vector<std::vector<double>> outData;
#if OUTER_NESTED_THREADS > 1
        #pragma omp parallel for num_threads(OUTER_NESTED_THREADS) default(none) shared(SAMPLES, coutMutex, curr, prgbar_segm, std::cout, start, end, step, N, H_List, S2_List, outData)
#endif
        for (int s = 1; s <= SAMPLES; s++) {
            std::vector<double> rawData = hlp::rungeKutta4_X(start, end, step, N, H_List, S2_List);

            coutMutex.lock();
            outData.emplace_back(rawData);
            int p = (int) ( (float) curr / (float) SAMPLES * (float) prgbar_segm);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < prgbar_segm; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) curr / (float) SAMPLES * 100.0 ) << "% (" << curr << "/" << SAMPLES << ")     ";
            std::cout.flush();
            curr++;
            coutMutex.unlock();
        }

        outData.shrink_to_fit();

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";
        std::cout << "preparing and saving data\n";

        ///// avg and saving results (for different sample sizes n) /////

        // gather x-data
        std::vector<double> beta_Data;
        double beta = start - step;
        while (beta <= end) {
            beta += step;
            beta_Data.emplace_back(beta);
        } beta_Data.shrink_to_fit();

        hlp::saveStatisticsData("data_susceptibility_J_const_QT",
                                "N: " + std::to_string(N) + "\n"
                                + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                                + "samples: " + std::to_string(SAMPLES) + "\n"
                                + "this took: " + formatTime(elapsed_seconds),
                                "beta in kb / J2", "C in J2", beta_Data, outData, SAMPLES, step, N);

    }

    ///// C and X /////

    void start_calculation_CX_J_const(const double &start, const double &end, const double &step, const double &J1,
                                      const double &J2, const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();
        std::cout << "\n" << "C(T) & X(T), J = const, QT, momentum states, N: " << N << ", step size: " << step << " ..." << std::endl;

        ///// get H and S2 /////
        std::vector<matrixTypeComplex> H_List;
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
            for (int k = k_lower; k <= k_upper; k++) {
                if (states.at(m).at(k - k_lower).empty()) {continue;}
                Eigen::MatrixXcd H_Mtrx = fillHamiltonBlock(J1, J2, 0.0, k, states.at(m).at(k - k_lower),
                                                            R_vals.at(m).at(k - k_lower), N);
                Eigen::MatrixXcd S2_Mtrx = ED::spinMatrixMomentum(N, k, states.at(m).at(k - k_lower), R_vals.at(m).at(k - k_lower));
                H_List.emplace_back(H_Mtrx.sparseView());
                S2_List.emplace_back(S2_Mtrx.sparseView());
            }
        }

        // init progressbar
        int prgbar_segm = 50;
        int curr = 0;
        coutMutex.lock();
        int _p = (int) ( (float) curr / (float) SAMPLES * (float) prgbar_segm);
        std::cout << "\r[";
        for (int _ = 0; _ < _p; _++) {
            std::cout << "#";
        } for (int _ = _p; _ < prgbar_segm; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) curr / (float) SAMPLES * 100.0 ) << "% (" << curr << "/" << SAMPLES << ")     ";
        std::cout.flush();
        curr++;
        coutMutex.unlock();

        std::vector<std::vector<double>> outDataC;
        std::vector<std::vector<double>> outDataX;


//#if OUTER_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(OUTER_NESTED_THREADS) default(none) shared(SAMPLES, coutMutex, curr, prgbar_segm, std::cout, start, end, step, N, H_List, S2_List, outDataC, outDataX)
//#endif
#pragma omp parallel for default(none) shared(SAMPLES, coutMutex, curr, prgbar_segm, std::cout, start, end, step, N, H_List, S2_List, outDataC, outDataX)
        for (int s = 1; s <= SAMPLES; s++) {
            std::vector<std::tuple<double, double>> rawData = hlp::rungeKutta4_CX(start, end, step, N, H_List, S2_List);
            std::vector<double> rawDataC;
            std::vector<double> rawDataX;
            for (std::tuple<double, double> tup : rawData) {
                rawDataC.emplace_back(std::get<0>(tup));
                rawDataX.emplace_back(std::get<1>(tup));
            }

            coutMutex.lock();
            outDataC.emplace_back(rawDataC);
            outDataX.emplace_back(rawDataX);
            int p = (int) ( (float) curr / (float) SAMPLES * (float) prgbar_segm);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < prgbar_segm; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) curr / (float) SAMPLES * 100.0 ) << "% (" << curr << "/" << SAMPLES << ")     ";
            std::cout.flush();
            curr++;
            coutMutex.unlock();
        }

        outDataC.shrink_to_fit();
        outDataX.shrink_to_fit();

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";
        std::cout << "preparing and saving data\n";

        ///// avg and saving results (for different sample sizes n) /////

        // gather x-data
        std::vector<double> beta_Data;
        double beta = start - step;
        while (beta <= end) {
            beta += step;
            beta_Data.emplace_back(beta);
        } beta_Data.shrink_to_fit();

        std::cout << "specific heat C:\n";
        hlp::saveStatisticsData("data_specific_heat_J_const_QT",
                                "N: " + std::to_string(N) + "\n"
                                + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                                + "h: " + std::to_string(0.0)+ "\n"
                                + "samples: " + std::to_string(SAMPLES) + "\n"
                                + "this took: " + formatTime(elapsed_seconds),
                                "beta in kb / J2", "C in J2", beta_Data, outDataC, SAMPLES, step, N);

        std::cout << "magnetic susceptibility X:\n";
        hlp::saveStatisticsData("data_susceptibility_J_const_QT",
                                "N: " + std::to_string(N) + "\n"
                                + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                                + "samples: " + std::to_string(SAMPLES) + "\n"
                                + "this took: " + formatTime(elapsed_seconds),
                                "beta in kb / J2", "C in J2", beta_Data, outDataX, SAMPLES, step, N);

    }

    ///// spin gap /////

    void start_calc_spin_gap(const double &J_START, const double &J_END, const int &J_COUNT,
                             const double &BETA_START, const double &BETA_END, const double &BETA_STEP,
                             const int &N, const int &SIZE, const int &SAMPLES) {

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
            for (int k = k_lower; k <= k_upper; k++) {
                if (states.at(m).at(k - k_lower).empty()) {continue;}
                Eigen::MatrixXcd S2_Mtrx = ED::spinMatrixMomentum(N, k, states.at(m).at(k - k_lower), R_vals.at(m).at(k - k_lower));
                S2_List.emplace_back(S2_Mtrx.sparseView());
            }
        }

        ///// gather x-data /////
        std::vector<double> beta_Data;
        double beta = BETA_START - BETA_STEP;
        while (beta <= BETA_END) {
            beta += BETA_STEP;
            beta_Data.emplace_back(beta);
        } beta_Data.shrink_to_fit();

        // init progressbar
        int prgbar_segm = 50;
        int curr = 0;
        coutMutex.lock();
        int _p = (int) ( (float) curr / (float) J_COUNT * (float) prgbar_segm);
        std::cout << "\r[";
        for (int _ = 0; _ < _p; _++) {
            std::cout << "#";
        } for (int _ = _p; _ < prgbar_segm; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) curr / (float) J_COUNT * 100.0 ) << "% (" << curr << "/" << J_COUNT << "), J1/J2 = " << J_START + (J_END - J_START) * curr / J_COUNT << "     ";
        std::cout.flush();
        curr++;
        coutMutex.unlock();

//#if OUTERMOST_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(OUTERMOST_NESTED_THREADS) default(none) shared(J_COUNT, J_START, J_END, N, SIZE, SAMPLES, coutMutex, BETA_START, BETA_END, BETA_STEP, S2_List, beta_Data, curr, prgbar_segm, std::cout)
//#endif
#pragma omp parallel for default(none) shared(J_COUNT, J_START, J_END, N, SIZE, SAMPLES, coutMutex, BETA_START, BETA_END, BETA_STEP, S2_List, beta_Data, curr, prgbar_segm, std::cout)
        for (int J_pos = 0; J_pos < J_COUNT; J_pos++) {
            double J = J_START + (J_END - J_START) * J_pos / J_COUNT;
            std::vector<matrixTypeComplex> H_List = getHamilton(J, 1.0, 0.0, N, SIZE);
            std::vector<std::vector<double>> rawDataX;
//#if OUTER_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(OUTER_NESTED_THREADS) default(none) shared(SAMPLES, coutMutex, BETA_START, BETA_END, BETA_STEP, S2_List, H_List, rawDataX, beta_Data, N, SIZE)
//#endif
            for (int s = 1; s <= SAMPLES; s++) {
                std::vector<double> rawData = hlp::rungeKutta4_X(BETA_START, BETA_END, BETA_STEP, N, H_List, S2_List);
                rawDataX.emplace_back(rawData);
            }
            // save data (silent)
            rawDataX.shrink_to_fit();
            hlp::saveAvgData("./results/" + std::to_string(N) + "/data/spin_gap_data/X_J" + std::to_string(J) + "QT.txt",
                             "samples: " + std::to_string(SAMPLES) + "\n",
                             "beta in kb / J2", "X in J2", beta_Data, rawDataX, N);
            // progressbar
            coutMutex.lock();
            int p = (int) ( (float) curr / (float) J_COUNT * (float) prgbar_segm);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < prgbar_segm; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) curr / (float) J_COUNT * 100.0 ) << "% (" << curr << "/" << J_COUNT << "), J1/J2 = " << J_START + (J_END - J_START) * curr / J_COUNT << "     ";
            std::cout.flush();
            curr++;
            coutMutex.unlock();
        }

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

    }

    ///// excitation energies /////

    void start_calc_excitation_energies(const double &J_START, const double &J_END, const int &J_COUNT,
                                        const double &BETA_START, const double &BETA_END, const double &BETA_STEP,
                                        const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();
        std::cout << "\n" << "excitation energies, QT, momentum states, N: " << N << ", step size: " << BETA_STEP << " ..." << std::endl;

        ///// gather x-data /////
        std::vector<double> beta_Data;
        double beta = BETA_START - BETA_STEP;
        while (beta <= BETA_END) {
            beta += BETA_STEP;
            beta_Data.emplace_back(beta);
        } beta_Data.shrink_to_fit();

        // init progressbar
        int prgbar_segm = 50;
        int curr = 0;
        coutMutex.lock();
        int _p = (int) ( (float) curr / (float) J_COUNT * (float) prgbar_segm);
        std::cout << "\r[";
        for (int _ = 0; _ < _p; _++) {
            std::cout << "#";
        } for (int _ = _p; _ < prgbar_segm; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) curr / (float) J_COUNT * 100.0 ) << "% (" << curr << "/" << J_COUNT << "), J1/J2 = " << J_START + (J_END - J_START) * curr / J_COUNT << "     ";
        std::cout.flush();
        curr++;
        coutMutex.unlock();

#pragma omp parallel for default(none) shared(J_COUNT, J_START, J_END, N, SIZE, SAMPLES, coutMutex, BETA_START, BETA_END, BETA_STEP, beta_Data, curr, prgbar_segm, std::cout)
        for (int J_pos = 0; J_pos < J_COUNT; J_pos++) {
            double J = J_START + (J_END - J_START) * J_pos / J_COUNT;
            std::vector<matrixTypeComplex> H_List = getHamilton(J, 1.0, 0.0, N, SIZE);
            std::vector<std::vector<double>> rawDataC;

            for (int s = 1; s <= SAMPLES; s++) {
                std::vector<double> rawData = hlp::rungeKutta4_C(BETA_START, BETA_END, BETA_STEP, N, H_List);
                rawDataC.emplace_back(rawData);
            }
            // save data (silent)
            rawDataC.shrink_to_fit();
            hlp::saveAvgData("./results/" + std::to_string(N) + "/data/excitation_energies_data/C_J" + std::to_string(J) + "QT.txt",
                             "samples: " + std::to_string(SAMPLES) + "\n",
                             "beta in kb / J2", "C in J2", beta_Data, rawDataC, N);
            // progressbar
            coutMutex.lock();
            int p = (int) ( (float) curr / (float) J_COUNT * (float) prgbar_segm);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < prgbar_segm; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) curr / (float) J_COUNT * 100.0 ) << "% (" << curr << "/" << J_COUNT << "), J1/J2 = " << J_START + (J_END - J_START) * curr / J_COUNT << "     ";
            std::cout.flush();
            curr++;
            coutMutex.unlock();
        }

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

    }

}
