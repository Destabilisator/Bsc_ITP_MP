#include "si.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::SI {

    ///// hamilton /////

    std::vector<matrixType> getHamilton(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE) {

        std::vector<matrixType> matList;

        std::vector<std::vector<int>> states_m(N+1);
        for (int s = 0; s < SIZE; s++) {
            states_m.at(ED::bitSum(s, N)).push_back(s);
        }

        std::vector<int> states;
        std::vector<int> R_vals;
        std::vector<int> m_vals;
        std::vector<int> n_vals;
        std::vector<int> c_vals;

        int numberOfStates = 0;
        const int k_upper = N/4;

        for (int mag = 0; mag <= N; mag++) {
            for (int k = 0; k <= k_upper; k++) {
                if (mag == N/2) {
                    for (int z : {-1, 1}) {
                        for (int p: {-1, 1}) {
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
                                Eigen::MatrixXd Mtrx = fillHamiltonSIBlock(J1, J2, h, k, p, z, states, R_vals, m_vals, n_vals, c_vals, N);
                                matList.emplace_back(Mtrx.sparseView());
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
                                    states.push_back(s);
                                    R_vals.push_back(sigma * R);
                                    m_vals.push_back(m);
                                    numberOfStates++;
                                }
                            }
                        }
                        if (!states.empty()) {
                            Eigen::MatrixXd Mtrx = fillHamiltonPSBlock(J1, J2, h, k, p, states, R_vals, m_vals, N);
                            matList.emplace_back(Mtrx.sparseView());
                        }
                        states.clear();
                        R_vals.clear();
                        m_vals.clear();
                    }
                }
            }
        }

        matList.shrink_to_fit();
        return matList;

    }

    Eigen::MatrixXd fillHamiltonPSBlock(const double &J1, const double &J2, const double &h, const int &k, const int &p,
                                        const std::vector<int> &states, const std::vector<int> &R_vals,
                                        const std::vector<int> &m_vals, const int &N) {

        int statesCount = (int) states.size();
        Eigen::MatrixXd hamiltonBlock = Eigen::MatrixXd::Zero(statesCount,statesCount);

        for (int a; a < statesCount; a++) {

            int s = states.at(a);

//            std::cout << hamiltonBlock(a,a) << std::endl;
            int mag = ED::bitSum(s, N) - N/2;
            hamiltonBlock(a,a) += (double) mag * h;
//            std::cout << hamiltonBlock(a,a) << std::endl << std::endl;

            int state_n = 1;
            if (a > 0 && states.at(a - 1) == states.at(a)) {
                continue;
            } else if (a < states.size() - 1 && states.at(a) == states.at(a + 1)) {
                state_n = 2;
            }

            for (int n = 0; n < N/2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J1;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -=  0.25 * J1;
                    }
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int r = 0, l = 0, q = 0;
                    ED::representative(d, &r, &l, &q, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J1 * ED::parityStates::helement(i, j, l, q, k, p, R_vals, m_vals, N);
                            }
                        }
                    }
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J2;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -= 0.25 * J2;
                    }
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    int r = 0, l = 0, q = 0;
                    ED::representative(d, &r, &l, &q, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J2 * ED::parityStates::helement(i, j, l, q, k, p, R_vals, m_vals, N);
                            }
                        }
                    }
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J2;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -= 0.25 * J2;
                    }
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    int r = 0, l = 0, q = 0;
                    ED::representative(d, &r, &l, &q, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J2 * ED::parityStates::helement(i, j, l, q, k, p, R_vals, m_vals, N);
                            }
                        }
                    }
                }
            }
        }

        return hamiltonBlock;

    }

    Eigen::MatrixXd fillHamiltonSIBlock(const double &J1, const double &J2, const double &h, const int &k, const int &p,
                                        const int &z, const std::vector<int> &states, const std::vector<int> &R_vals,
                                        const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                                        const std::vector<int> &c_vals, const int &N) {

        int statesCount = (int) states.size();
        Eigen::MatrixXd hamiltonBlock = Eigen::MatrixXd::Zero(statesCount,statesCount);

        for (int a; a < statesCount; a++) {

            int s = states.at(a);

            // external magnetic field
//            std::cout << hamiltonBlock(a,a) << std::endl;
//            int mag = bitSum(s, N) - N/2;
//            hamiltonBlock(a,a) += (double) mag * h;
//            std::cout << hamiltonBlock(a,a) << std::endl << std::endl;

            int state_n = 1;
            if (a > 0 && states.at(a - 1) == states.at(a)) {
                continue;
            } else if (a < states.size() - 1 && states.at(a) == states.at(a + 1)) {
                state_n = 2;
            }

            for (int n = 0; n < N/2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J1;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -=  0.25 * J1;
                    }
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int r = 0, l = 0, q = 0, g = 0;
                    ED::representative(d, r, l, q, g, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J1 * ED::spinInversion::helement(i, j, l, q, g, k, p, z, R_vals, m_vals, n_vals, c_vals, N);
                            }
                        }
                    }
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J2;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -= 0.25 * J2;
                    }
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    int r = 0, l = 0, q = 0, g = 0;
                    ED::representative(d, r, l, q, g, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J2 * ED::spinInversion::helement(i, j, l, q, g, k, p, z, R_vals, m_vals, n_vals, c_vals, N);
                            }
                        }
                    }
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J2;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -= 0.25 * J2;
                    }
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    int r = 0, l = 0, q = 0, g = 0;
                    ED::representative(d, r, l, q, g, N);
                    int b = ED::findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J2 * ED::spinInversion::helement(i, j, l, q, g, k, p, z, R_vals, m_vals, n_vals, c_vals, N);
                            }
                        }
                    }
                }
            }
        }

        return hamiltonBlock;

    }


    ///// C /////
/*
    void start_calculation_C_J_const(const double &start, const double &end, const double &step, const double &J1,
                                     const double &J2, const double &h, const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();
        std::cout << "\n" << "C(T), J = const, QT, momentum states, N: " << N << ", step size: " << step << " ..." << std::endl;

        std::vector<matrixType> matrixList = getHamilton(J1, J2, h, N, SIZE);

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
/*
    ///// X /////

    void start_calculation_X_J_const(const double &start, const double &end, const double &step, const double &J1,
                                     const double &J2, const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();
        std::cout << "\n" << "X(T), J = const, QT, momentum states, N: " << N << ", step size: " << step << " ..." << std::endl;

        ///// get H and S2 /////
        std::vector<matrixType> H_List;
        std::vector<matrixType> S2_List;
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
        std::vector<matrixType> H_List;
        std::vector<matrixType> S2_List;
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

#if OUTER_NESTED_THREADS > 1
#pragma omp parallel for num_threads(OUTER_NESTED_THREADS) default(none) shared(SAMPLES, coutMutex, curr, prgbar_segm, std::cout, start, end, step, N, H_List, S2_List, outDataC, outDataX)
#endif
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
        std::vector<matrixType> S2_List;
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

#if OUTERMOST_NESTED_THREADS > 1
#pragma omp parallel for num_threads(OUTERMOST_NESTED_THREADS) default(none) shared(J_COUNT, J_START, J_END, N, SIZE, SAMPLES, coutMutex, BETA_START, BETA_END, BETA_STEP, S2_List, beta_Data, curr, prgbar_segm, std::cout)
#endif
        for (int J_pos = 0; J_pos < J_COUNT; J_pos++) {
            double J = J_START + (J_END - J_START) * J_pos / J_COUNT;
            std::vector<matrixType> H_List = getHamilton(J, 1.0, 0.0, N, SIZE);
            std::vector<std::vector<double>> rawDataX;
#if OUTER_NESTED_THREADS > 1
#pragma omp parallel for num_threads(OUTER_NESTED_THREADS) default(none) shared(SAMPLES, coutMutex, BETA_START, BETA_END, BETA_STEP, S2_List, H_List, rawDataX, beta_Data, N, SIZE)
#endif
            for (int s = 1; s <= SAMPLES; s++) {
                std::vector<double> rawData = hlp::rungeKutta4_X(BETA_START, BETA_END, BETA_STEP, N, H_List, S2_List);
                rawDataX.emplace_back(rawData);
            }
            // save data (silent)
            rawDataX.shrink_to_fit();
            hlp::saveAvgData("./results/" + std::to_string(N) + "/data/spin_gap_data/X_J" + std::to_string(J) + ".txt",
                             "samples: " + std::to_string(SAMPLES) + "\n",
                             "beta in kb / J2", "C in J2", beta_Data, rawDataX, N);
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
*/
}
