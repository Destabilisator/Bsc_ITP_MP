#include "hamilton.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace QT::hamilton {

    Eigen::SparseMatrix<std::complex<double>> getHamilton(const double &J1, const double &J2, const int &N, const int &SIZE) {

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

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                fillHamiltonBlockMomentum(J1, J2, k, states.at(m).at(k - k_lower), R_vals.at(m).at(k - k_lower),
                                  matrixBlocks, N, SIZE);
            }
        }

        int offset_blocks = 0;
        // fill Triplets
        typedef Eigen::Triplet<std::complex<double>> Trip;
        std::vector<Trip> trp;
        for (const Eigen::MatrixXcd &M: matrixBlocks) {
            for (int i = 0; i < (int) M.rows(); i++) {
                for (int j = 0; j < (int) M.rows(); j++) {
                    trp.emplace_back(offset_blocks + i, offset_blocks + j, M(i,j));
                }
            }
            offset_blocks += (int) M.rows();
        }

        // create SparseMatrix from Triplets
        Eigen::SparseMatrix<std::complex<double>> H(SIZE, SIZE);
        H.setFromTriplets(trp.begin(), trp.end());
//        H.makeCompressed();
        return H;
    }

    void fillHamiltonBlockMomentum(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                           const std::vector<int> &R_vals, std::vector<Eigen::MatrixXcd> &hamiltonBlocks,
                           const int &N, const int &SIZE) {

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
        hamiltonBlocks.emplace_back(hamiltonBlock);
    }

}
