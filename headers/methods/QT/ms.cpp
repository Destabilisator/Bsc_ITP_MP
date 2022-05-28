#include "ms.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::MS {

    std::vector<matrixDataMomentumType> getHamilton(const double &J1, const double &J2, const int &N, const int &SIZE) {

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
                                                          R_vals.at(m).at(k - k_lower), N, SIZE);
                matList.emplace_back(m, k, Mtrx.sparseView());
            }
        }

        matList.shrink_to_fit();
        return matList;

    }

    Eigen::MatrixXcd fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                                       const std::vector<int> &R_vals, const int &N, const int &SIZE) {

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

    std::vector<Eigen::VectorXcd> getVector(const int &N, const int &SIZE, const std::vector<matrixDataMomentumType> &matrixBlocks) {

        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        std::vector<Eigen::VectorXcd> vectors;

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

}
