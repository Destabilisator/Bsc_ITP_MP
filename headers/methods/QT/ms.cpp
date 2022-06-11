#include "ms.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::MS {

    Eigen::VectorXcd getVector(int size) {

        // random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        // fill vector
        Eigen::VectorXcd v = Eigen::VectorXcd::Zero(size);
        for (int i = 0; i < size; i++) {
            v(i) = std::complex<double>(randNum(gen), randNum(gen));
        }
        v.normalize();
        return v;

    }

    std::vector<Eigen::VectorXcd> getVector(const int &N, const int &SIZE, const std::vector<matrixDataMomentumType> &matrixBlocks) {

        // random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        std::vector<Eigen::VectorXcd> vectors;

        // fill vector
        for (const matrixDataMomentumType &data : matrixBlocks) {
            int size = (int) std::get<2>(data).size();
            Eigen::VectorXcd v = Eigen::VectorXcd::Zero(size);
            for (int i = 0; i < size; i++) {
                v(i) = std::complex<double>(randNum(gen), randNum(gen));
            }
            v.normalize();
            vectors.push_back(v);
            std::cout << v << std::endl;
        }

        return vectors;

    }

    std::vector<Eigen::VectorXcd> getVector(const std::vector<matrixType> &matrixBlocks) {

        // random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        std::vector<Eigen::VectorXcd> vectors;

        // fill vector
        for (const matrixType &block : matrixBlocks) {
            int size = (int) block.rows();
            Eigen::VectorXcd v = Eigen::VectorXcd::Zero(size);
            for (int i = 0; i < size; i++) {
                v(i) = std::complex<double>(randNum(gen), randNum(gen));
            }
            v.normalize();
            vectors.push_back(v);
//            std::cout << v << std::endl;
        }

        return vectors;

    }

    ///// C /////

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

    std::vector<matrixType> getHamilton(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE) {

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

        std::vector<matrixType> matList;

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

    std::vector<double> rungeKutta4_C(const double &start, const double &end, const double &step, const int &N,
                                      const std::vector<matrixType> &matrixList) {

        std::vector<double> outData;
        std::vector<Eigen::VectorXcd> vec = getVector(matrixList);
        int blockCount = (int) matrixList.size();

        double norm = 0.0;
        for (int i = 0; i < blockCount; i++) {
            norm += std::pow(vec.at(i).norm(), 2);
        } norm = std::sqrt(norm);
        for (int i = 0; i < blockCount; i++) {
            vec.at(i) = vec.at(i) / norm;
        }

        double beta = start - step;
        while (beta <= end) {

            beta += step;
            norm = 0.0;
#if INNER_NESTED_THREADS > 1
            #pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, vec, matrixList, step, norm)
#endif
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = hlp::rungeKutta4Block(vec.at(i), matrixList.at(i), step);
                double normnt = std::pow(vec.at(i).norm(), 2);
                #pragma omp critical
                norm += normnt;
            }

            norm = std::sqrt(norm);
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = vec.at(i) / norm;
            }

            std::vector<double> vec_H_vec_List;
            std::vector<double> vec_H2_vec_List;

            //double C = 0.0;
#if INNER_NESTED_THREADS > 1
            #pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, matrixList, vec, vec_H_vec_List, vec_H2_vec_List)
#endif
            for (int i = 0; i < blockCount; i++) {
                const matrixType &H = matrixList.at(i);
                const Eigen::VectorXcd &v = vec.at(i);
                double vHv = std::real((v.adjoint() * H * v)(0,0));
                double vH2v = std::real((v.adjoint() * H * H * v)(0,0));
                #pragma omp critical
                vec_H_vec_List.emplace_back(vHv);
                #pragma omp critical
                vec_H2_vec_List.emplace_back(vH2v);
            }

            double vec_H_vec = std::accumulate(vec_H_vec_List.begin(), vec_H_vec_List.end(), 0.0);
            double vec_H2_vec = std::accumulate(vec_H2_vec_List.begin(), vec_H2_vec_List.end(), 0.0);

            double H_diff = std::real(vec_H2_vec - std::pow(vec_H_vec, 2));
            double C = beta * beta * H_diff / (double) N;

            outData.emplace_back(C);

        }

        outData.shrink_to_fit();
        return outData;

    }

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
            std::vector<double> rawData = rungeKutta4_C(start, end, step, N, matrixList);
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

/*
#ifdef SAVE_WITH_SETS_OF_n_SAMPLES
        // save for different amounts n of random states
        for (int n = 1; n <= SAMPLES / 2; n++) {
            std::vector<std::vector<double>> C_Data_to_avg;
            // combine data for every i to i+n samples into one avg
            for (int i = 0; i+n < SAMPLES; i += n) {
                std::vector<std::vector<double>> C_Data_raw;
                // get data to avg
                for (int j = i; j < i + n; j++) {
//                    std:: cout << "n: " << n << ", i: " << i << ", j: " << j << std::endl;
                    C_Data_raw.emplace_back(outData.at(j));
                }
                std::vector<double> C_vector;
                // avg
                for(int k = 0; k < C_Data_raw.at(0).size(); k++) {
                    std::vector<double> C_vector_raw;
                    for (std::vector<double> Cd : C_Data_raw) {
                        C_vector_raw.emplace_back(Cd.at(k));
                    }
                    std::tuple<double, double> mean_se = get_mean_and_se(C_vector_raw);
                    C_vector.emplace_back(std::get<0>(mean_se));
                }
                C_Data_to_avg.emplace_back(C_vector);
            }
            std::vector<double> C_Data_save;
            std::vector<double> CErr_Data_save;
            // avg and stdv of avg C with n samples each
            for(int k = 0; k < C_Data_to_avg.at(0).size(); k++) {
                std::vector<double> C_vector_raw;
                for (std::vector<double> Cd : C_Data_to_avg) {
                    C_vector_raw.emplace_back(Cd.at(k));
                } std::tuple<double, double> mean_se_C = get_mean_and_se(C_vector_raw);
                C_Data_save.emplace_back(std::get<0>(mean_se_C));
                CErr_Data_save.emplace_back(std::get<1>(mean_se_C));
            }
            // save data N_n_data_specific_heat_J_const_QT_txt
            std::string filename = std::to_string(n) + "_data_specific_heat_J_const_QT";
            hlp::saveOutData(filename + ".txt",
                             "N: " + std::to_string(N) + "\n"
                             + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                             + "h: " + std::to_string(h)+ "\n"
                             + "samples: " + std::to_string(SAMPLES) + "\n"
                             + "this took: " + formatTime(elapsed_seconds),
                             "beta in kb / J2", "C in J2", beta_Data, C_Data_save, CErr_Data_save, N);
#ifdef SAVE_WITH_STEP_SIZE
            hlp::saveOutData("step_size_data/" + filename + "_step" + std::to_string(step) + ".txt",
                             "N: " + std::to_string(N) + "\n"
                             + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                             + "h: " + std::to_string(h)+ "\n"
                             + "samples: " + std::to_string(SAMPLES) + "\n"
                             + "this took: " + formatTime(elapsed_seconds),
                             "beta in kb / J2", "C in J2", beta_Data, C_Data_save, CErr_Data_save, N);
#endif

        }

#endif

#ifdef SAVE_WITH_DATA_FROM_ALL_SAMPLES
        std::vector<double> C_Data;
        std::vector<double> CErr_Data;

        // avg and stdv
        for(int i = 0; i < outData.at(0).size(); i++) {
            std::vector<double> C_temp_data;
            for (std::vector<double> C_data_raw : outData) {
                C_temp_data.emplace_back(C_data_raw.at(i));
            }
            std::tuple<double, double> mean_se = get_mean_and_se(C_temp_data);
            C_Data.emplace_back(std::get<0>(mean_se));
            CErr_Data.emplace_back(std::get<1>(mean_se));
        }
        std::string filename = "data_specific_heat_J_const_QT";
        hlp::saveOutData(filename + ".txt",
                         "N: " + std::to_string(N) + "\n"
                         + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                         + "h: " + std::to_string(h)+ "\n"
                         + "samples: " + std::to_string(SAMPLES) + "\n"
                         + "this took: " + formatTime(elapsed_seconds),
                         "beta in kb / J2", "C in J2", beta_Data, C_Data, CErr_Data, N);
#ifdef SAVE_WITH_STEP_SIZE
        hlp::saveOutData("step_size_data/" + filename + "_step" + std::to_string(step) + ".txt",
                         "N: " + std::to_string(N) + "\n"
                         + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                         + "h: " + std::to_string(h)+ "\n"
                         + "samples: " + std::to_string(SAMPLES) + "\n"
                         + "this took: " + formatTime(elapsed_seconds),
                         "beta in kb / J2", "C in J2", beta_Data, C_Data, CErr_Data, N);
#endif

#endif
*/
    }

    ///// X /////

    std::vector<double> rungeKutta4_X(const double &start, const double &end, const double &step, const int &N,
                                      const std::vector<matrixType> &H_List, const std::vector<matrixType> &S2_List) {

        std::vector<double> outData;
        std::vector<Eigen::VectorXcd> vec = getVector(H_List);
        int blockCount = (int) H_List.size();

        double norm = 0.0;
        for (int i = 0; i < blockCount; i++) {
            norm += std::pow(vec.at(i).norm(), 2);
        } norm = std::sqrt(norm);
        for (int i = 0; i < blockCount; i++) {
            vec.at(i) = vec.at(i) / norm;
        }

        double beta = start - step;
        while (beta <= end) {

            beta += step;
            norm = 0.0;

            // RK4 with vec on H to get new state
#if INNER_NESTED_THREADS > 1
#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, vec, H_List, step, norm)
#endif
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = hlp::rungeKutta4Block(vec.at(i), H_List.at(i), step);
                double normnt = std::pow(vec.at(i).norm(), 2);
#pragma omp critical
                norm += normnt;
            }

            // ensure norm(vec) = 1
            norm = std::sqrt(norm);
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = vec.at(i) / norm;
            }

            std::vector<double> vec_S2_vec_List;

#if INNER_NESTED_THREADS > 1
#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, S2_List, vec, vec_S2_vec_List)
#endif
            for (int i = 0; i < blockCount; i++) {
                const matrixType &S2 = S2_List.at(i);
                const Eigen::VectorXcd &v = vec.at(i);
                double vS2v = std::real((v.adjoint() * S2 * v)(0,0));
#pragma omp critical
                vec_S2_vec_List.emplace_back(vS2v);
            }

            double vec_S2_vec = std::accumulate(vec_S2_vec_List.begin(), vec_S2_vec_List.end(), 0.0);

            double X = beta * vec_S2_vec / 3.0 / (double) N;

            outData.emplace_back(X);

        }

        outData.shrink_to_fit();
        return outData;

    }

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
            std::vector<double> rawData = rungeKutta4_X(start, end, step, N, H_List, S2_List);

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

/*
#ifdef SAVE_WITH_SETS_OF_n_SAMPLES
        // save for different amounts n of random states
        for (int n = 1; n <= SAMPLES / 2; n++) {
            std::vector<std::vector<double>> X_Data_to_avg;
            // combine data for every i to i+n samples into one avg
            for (int i = 0; i+n < SAMPLES; i += n) {
                std::vector<std::vector<double>> X_Data_raw;
                // get data to avg
                for (int j = i; j < i + n; j++) {
//                    std:: cout << "n: " << n << ", i: " << i << ", j: " << j << std::endl;
                    X_Data_raw.emplace_back(outData.at(j));
                }
                std::vector<double> X_vector;
                // avg
                for(int k = 0; k < X_Data_raw.at(0).size(); k++) {
                    std::vector<double> X_vector_raw;
                    for (std::vector<double> Xd : X_Data_raw) {
                        X_vector_raw.emplace_back(Xd.at(k));
                    }
                    std::tuple<double, double> mean_se = get_mean_and_se(X_vector_raw);
                    X_vector.emplace_back(std::get<0>(mean_se));
                }
                X_Data_to_avg.emplace_back(X_vector);
            }
            std::vector<double> X_Data_save;
            std::vector<double> XErr_Data_save;
            // avg and stdv of avg C with n samples each
            for(int k = 0; k < X_Data_to_avg.at(0).size(); k++) {
                std::vector<double> X_vector_raw;
                for (std::vector<double> Xd : X_Data_to_avg) {
                    X_vector_raw.emplace_back(Xd.at(k));
                } std::tuple<double, double> mean_se_X = get_mean_and_se(X_vector_raw);
                X_Data_save.emplace_back(std::get<0>(mean_se_X));
                XErr_Data_save.emplace_back(std::get<1>(mean_se_X));
            }
            // save data N_n_data_specific_heat_J_const_QT_txt
            std::string filename = std::to_string(n) + "_data_susceptibility_J_const_QT";
            hlp::saveOutData(filename + ".txt",
                             "N: " + std::to_string(N) + "\n"
                             + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                             + "samples: " + std::to_string(SAMPLES) + "\n"
                             + "this took: " + formatTime(elapsed_seconds),
                             "beta in kb / J2", "C in J2", beta_Data, X_Data_save, XErr_Data_save, N);
#ifdef SAVE_WITH_STEP_SIZE
            hlp::saveOutData("step_size_data/" + filename + "_step" + std::to_string(step) + ".txt",
                             "N: " + std::to_string(N) + "\n"
                             + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                             + "samples: " + std::to_string(SAMPLES) + "\n"
                             + "this took: " + formatTime(elapsed_seconds),
                             "beta in kb / J2", "C in J2", beta_Data, X_Data_save, XErr_Data_save, N);
#endif

        }

#endif

#ifdef SAVE_WITH_DATA_FROM_ALL_SAMPLES
        std::vector<double> X_Data;
        std::vector<double> XErr_Data;

        // avg and stdv of X
        for(int i = 0; i < outData.at(0).size(); i++) {
            std::vector<double> X_temp_data;
            for (std::vector<double> C_data_raw : outData) {
                X_temp_data.emplace_back(C_data_raw.at(i));
            }
            std::tuple<double, double> mean_se = get_mean_and_se(X_temp_data);
            X_Data.emplace_back(std::get<0>(mean_se));
            XErr_Data.emplace_back(std::get<1>(mean_se));
        }
        std::string filename = "data_susceptibility_J_const_QT";
        hlp::saveOutData(filename + ".txt",
                         "N: " + std::to_string(N) + "\n"
                         + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                         + "samples: " + std::to_string(SAMPLES) + "\n"
                         + "this took: " + formatTime(elapsed_seconds),
                         "beta in kb / J2", "C in J2", beta_Data, X_Data, XErr_Data, N);
#ifdef SAVE_WITH_STEP_SIZE
        hlp::saveOutData("step_size_data/" + filename + "_step" + std::to_string(step) + ".txt",
                         "N: " + std::to_string(N) + "\n"
                         + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                         + "samples: " + std::to_string(SAMPLES) + "\n"
                         + "this took: " + formatTime(elapsed_seconds),
                         "beta in kb / J2", "C in J2", beta_Data, X_Data, XErr_Data, N);
#endif

#endif
*/
    }

    ///// C_X /////

    std::vector<std::tuple<double, double>> rungeKutta4_CX(const double &start, const double &end, const double &step, const int &N,
                                      const std::vector<matrixType> &H_List, const std::vector<matrixType> &S2_List) {

        std::vector<std::tuple<double, double>> outData;
        std::vector<Eigen::VectorXcd> vec = getVector(H_List);
        int blockCount = (int) H_List.size();

        double norm = 0.0;
        for (int i = 0; i < blockCount; i++) {
            norm += std::pow(vec.at(i).norm(), 2);
        } norm = std::sqrt(norm);
        for (int i = 0; i < blockCount; i++) {
            vec.at(i) = vec.at(i) / norm;
        }

        double beta = start - step;
        while (beta <= end) {

            beta += step;
            norm = 0.0;

            // RK4 with vec on H to get new state
#if INNER_NESTED_THREADS > 1
#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, vec, H_List, step, norm)
#endif
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = hlp::rungeKutta4Block(vec.at(i), H_List.at(i), step);
                double normnt = std::pow(vec.at(i).norm(), 2);
#pragma omp critical
                norm += normnt;
            }

            // ensure norm(vec) = 1
            norm = std::sqrt(norm);
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = vec.at(i) / norm;
            }

            std::vector<double> vec_H_vec_List;
            std::vector<double> vec_H2_vec_List;
            std::vector<double> vec_S2_vec_List;

#if INNER_NESTED_THREADS > 1
#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, H_List, S2_List, vec, vec_S2_vec_List, vec_H_vec_List, vec_H2_vec_List)
#endif
            for (int i = 0; i < blockCount; i++) {
                const matrixType &H = H_List.at(i);
                const matrixType &S2 = S2_List.at(i);
                const Eigen::VectorXcd &v = vec.at(i);
                double vS2v = std::real((v.adjoint() * S2 * v)(0,0));
                double vHv = std::real((v.adjoint() * H * v)(0,0));
                double vH2v = std::real((v.adjoint() * H * H * v)(0,0));
#pragma omp critical
                vec_S2_vec_List.emplace_back(vS2v);
#pragma omp critical
                vec_H_vec_List.emplace_back(vHv);
#pragma omp critical
                vec_H2_vec_List.emplace_back(vH2v);
            }

            double vec_H_vec = std::accumulate(vec_H_vec_List.begin(), vec_H_vec_List.end(), 0.0);
            double vec_H2_vec = std::accumulate(vec_H2_vec_List.begin(), vec_H2_vec_List.end(), 0.0);
            double H_diff = std::real(vec_H2_vec - std::pow(vec_H_vec, 2));
            double C = beta * beta * H_diff / (double) N;

            double vec_S2_vec = std::accumulate(vec_S2_vec_List.begin(), vec_S2_vec_List.end(), 0.0);
            double X = beta * vec_S2_vec / 3.0 / (double) N;

            outData.emplace_back(C,X);

        }

        outData.shrink_to_fit();
        return outData;

    }

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
            std::vector<std::tuple<double, double>> rawData = rungeKutta4_CX(start, end, step, N, H_List, S2_List);
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

/*
#ifdef SAVE_WITH_SETS_OF_n_SAMPLES
        // save for different amounts n of random states
        for (int n = 1; n <= SAMPLES / 2; n++) {
            std::vector<std::vector<double>> X_Data_to_avg;
            // combine data for every i to i+n samples into one avg
            for (int i = 0; i+n < SAMPLES; i += n) {
                std::vector<std::vector<double>> X_Data_raw;
                // get data to avg
                for (int j = i; j < i + n; j++) {
//                    std:: cout << "n: " << n << ", i: " << i << ", j: " << j << std::endl;
                    X_Data_raw.emplace_back(outData.at(j));
                }
                std::vector<double> X_vector;
                // avg
                for(int k = 0; k < X_Data_raw.at(0).size(); k++) {
                    std::vector<double> X_vector_raw;
                    for (std::vector<double> Xd : X_Data_raw) {
                        X_vector_raw.emplace_back(Xd.at(k));
                    }
                    std::tuple<double, double> mean_se = get_mean_and_se(X_vector_raw);
                    X_vector.emplace_back(std::get<0>(mean_se));
                }
                X_Data_to_avg.emplace_back(X_vector);
            }
            std::vector<double> X_Data_save;
            std::vector<double> XErr_Data_save;
            // avg and stdv of avg C with n samples each
            for(int k = 0; k < X_Data_to_avg.at(0).size(); k++) {
                std::vector<double> X_vector_raw;
                for (std::vector<double> Xd : X_Data_to_avg) {
                    X_vector_raw.emplace_back(Xd.at(k));
                } std::tuple<double, double> mean_se_X = get_mean_and_se(X_vector_raw);
                X_Data_save.emplace_back(std::get<0>(mean_se_X));
                XErr_Data_save.emplace_back(std::get<1>(mean_se_X));
            }
            // save data N_n_data_specific_heat_J_const_QT_txt
            std::string filename = std::to_string(n) + "_data_susceptibility_J_const_QT";
            hlp::saveOutData(filename + ".txt",
                             "N: " + std::to_string(N) + "\n"
                             + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                             + "samples: " + std::to_string(SAMPLES) + "\n"
                             + "this took: " + formatTime(elapsed_seconds),
                             "beta in kb / J2", "C in J2", beta_Data, X_Data_save, XErr_Data_save, N);
#ifdef SAVE_WITH_STEP_SIZE
            hlp::saveOutData("step_size_data/" + filename + "_step" + std::to_string(step) + ".txt",
                             "N: " + std::to_string(N) + "\n"
                             + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                             + "samples: " + std::to_string(SAMPLES) + "\n"
                             + "this took: " + formatTime(elapsed_seconds),
                             "beta in kb / J2", "C in J2", beta_Data, X_Data_save, XErr_Data_save, N);
#endif

        }

#endif

#ifdef SAVE_WITH_DATA_FROM_ALL_SAMPLES
        std::vector<double> X_Data;
        std::vector<double> XErr_Data;

        // avg and stdv of X
        for(int i = 0; i < outData.at(0).size(); i++) {
            std::vector<double> X_temp_data;
            for (std::vector<double> C_data_raw : outData) {
                X_temp_data.emplace_back(C_data_raw.at(i));
            }
            std::tuple<double, double> mean_se = get_mean_and_se(X_temp_data);
            X_Data.emplace_back(std::get<0>(mean_se));
            XErr_Data.emplace_back(std::get<1>(mean_se));
        }
        std::string filename = "data_susceptibility_J_const_QT";
        hlp::saveOutData(filename + ".txt",
                         "N: " + std::to_string(N) + "\n"
                         + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                         + "samples: " + std::to_string(SAMPLES) + "\n"
                         + "this took: " + formatTime(elapsed_seconds),
                         "beta in kb / J2", "C in J2", beta_Data, X_Data, XErr_Data, N);
#ifdef SAVE_WITH_STEP_SIZE
        hlp::saveOutData("step_size_data/" + filename + "_step" + std::to_string(step) + ".txt",
                         "N: " + std::to_string(N) + "\n"
                         + "J1/J2: " + std::to_string(J1/J2)+ "\n"
                         + "samples: " + std::to_string(SAMPLES) + "\n"
                         + "this took: " + formatTime(elapsed_seconds),
                         "beta in kb / J2", "C in J2", beta_Data, X_Data, XErr_Data, N);
#endif

#endif
*/
    }

    ///// spin gap /////

    void start_calc_spin_gap(const double &J_START, const double &J_END, const int &J_COUNT,
                             const double &BETA_START, const double &BETA_END, const double &BETA_STEP,
                             const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();
        std::cout << "\n" << "C(T) J = const, QT, momentum states, N: " << N << ", step size: " << BETA_STEP << " ..." << std::endl;

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
                std::vector<double> rawData = rungeKutta4_X(BETA_START, BETA_END, BETA_STEP, N, H_List, S2_List);
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
}
