#include "ms.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::MS {

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

    std::vector<matrixDataMomentumType> getIndexAndHamilton(const double &J1, const double &J2, const int &N, const int &SIZE) {

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
                Eigen::MatrixXcd Mtrx = fillHamiltonBlock(J1, J2, k, states.at(m).at(k - k_lower),
                                                          R_vals.at(m).at(k - k_lower), N);
                matList.emplace_back(m, k, Mtrx.sparseView());
            }
        }

        matList.shrink_to_fit();
        return matList;

    }

    std::vector<matrixType> getHamilton(const double &J1, const double &J2, const int &N, const int &SIZE) {

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
        std::vector<matrixType> matList;

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                if (states.at(m).at(k - k_lower).empty()) {continue;}
                Eigen::MatrixXcd Mtrx = fillHamiltonBlock(J1, J2, k, states.at(m).at(k - k_lower),
                                                          R_vals.at(m).at(k - k_lower), N);
                matList.emplace_back(Mtrx.sparseView());
            }
        }

        matList.shrink_to_fit();
        return matList;

    }

    Eigen::MatrixXcd fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                                       const std::vector<int> &R_vals, const int &N) {

        const int statesCount = (int) states.size();
        Eigen::MatrixXcd hamiltonBlock = Eigen::MatrixXcd::Zero(statesCount,statesCount);

        for (int a = 0; a < statesCount; a++) {
            int s = states.at(a);
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

    std::vector<std::tuple<double, double>> rungeKutta4_C(const double &start, const double &end, const double &step,
                                                          const double &J1, const double &J2, const int &N, const int &SIZE,
                                                          const std::vector<matrixType> &matrixList) {

        std::vector<std::tuple<double, double>> outData;
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

            #pragma omp parallel for default(none) shared(blockCount, vec, matrixList, step, norm)
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = hlp::rungeKutta4Block(vec.at(i), matrixList.at(i), step);
                #pragma omp critical
                norm += std::pow(vec.at(i).norm(), 2);
            }

            norm = std::sqrt(norm);
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = vec.at(i) / norm;
            }

            std::vector<double> vec_H_vec_List;
            std::vector<double> vec_H2_vec_List;

            double C = 0.0;
            #pragma omp parallel for default(none) shared(blockCount, matrixList, vec, vec_H_vec_List, vec_H2_vec_List)
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
            C += beta * beta * H_diff / (double) N;

            outData.emplace_back(beta, C);

        }

        outData.shrink_to_fit();
        return outData;

    }

    void start_calculation_C_J_const(const double &start, const double &end, const double &step,
                            const double &J1, const double &J2, const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();

        std::cout << "\n" << "C(T), J = const, QT, momentum states ..." << std::endl;

        std::vector<matrixType> matrixList = getHamilton(J1, J2, N, SIZE);

        typedef std::vector<std::tuple<double, double>> dataVectorType;
        std::vector<dataVectorType> dataToAvg;
        for (int _ = 1; _ <= SAMPLES; _++) {
            std::cout << _ << "/" << SAMPLES << "\n";
            dataVectorType rawData = rungeKutta4_C(start, end, step, J1, J2, N, SIZE, matrixList);
            dataToAvg.push_back(rawData);
        }

        std::vector<std::tuple<double, double>> outData;

        for (int i = 0; i < dataToAvg.at(0).size(); i++) {
            double C = 0.0;
            double beta = std::get<0>(dataToAvg.at(0).at(i));
            for (dataVectorType data : dataToAvg) {
                C += std::get<1>(data.at(i));
            }
            outData.emplace_back(beta, C / (double) SAMPLES);
        }

        hlp::saveOutData("data_specific_heat_J_const_QT.txt", "QT, MS für N = " + std::to_string(N)
                        + " mit " + std::to_string(SAMPLES) + " Samples",
                        "T in J2 / kb", "C in J2", outData, N);

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

    }

}
