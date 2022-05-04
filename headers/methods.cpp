#include "methods.h"
#include "helpers.h"

/////////////////////////////// naiver Ansatz ///////////////////////////////
namespace naiv {
    void fillHamilton(double **hamilton, const double &J1,const double &J2, const int &N, const int &SIZE) {
        for (int s = 0; s <= SIZE - 1; s++) {
            for (int n = 0; n < N / 2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    hamilton[s][s] += 0.25 * J1;
                } else {
                    hamilton[s][s] -= 0.25 * J1;
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    hamilton[s][d] = 0.5 * J1;
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    hamilton[s][s] += 0.25 * J2;
                } else {
                    hamilton[s][s] -= 0.25 * J2;
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    hamilton[s][d] = 0.5 * J2;
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    hamilton[s][s] += 0.25 * J2;
                } else {
                    hamilton[s][s] -= 0.25 * J2;
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    hamilton[s][d] = 0.5 * J2;
                }
            }
        }
    }

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   const int &N, const int &SIZE) {
        static auto **hamilton1 = new double*[SIZE];
        for (int i = 0; i < SIZE; i++) {
            hamilton1[i] = new double[SIZE];
            for (int j = 0; j < SIZE; j++) {
                hamilton1[i][j] = 0.0;
            }
        }
        fillHamilton(hamilton1, J1, J2, N, SIZE);
        Eigen::MatrixXd H1(SIZE, SIZE);
        for (int i = 0; i < SIZE; i++) {
            H1.row(i) = Eigen::VectorXd::Map(&hamilton1[i][0], SIZE);
        }
#ifdef showMatrix
        std::cout << H1 << std::endl;
#endif
#ifdef saveMatrix
        saveHamilton(hamilton1, "HamiltonNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), SIZE);
#endif
        //std::cout << "solving...\n";
        Eigen::EigenSolver<Eigen::MatrixXd> solver1(H1);
        const Eigen::VectorXcd &H1EiVal = solver1.eigenvalues();
        for (std::complex<double> ev : H1EiVal) {
            HEiValList->push_back(ev);
        }
        // sort eigenvalues
#if defined(showEigenvalues) || defined(saveEigenvalues)
        HEiValList->shrink_to_fit();
        std::sort(HEiValList->begin(), HEiValList->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });
#endif
#ifdef showEigenvalues
        std::cout << "eigenvalues:\n";
        for (std::complex<double> ev : *HEiValList) {
            std::cout << ev << "\n";
        }
#endif
#ifdef saveEigenvalues
        saveComplexEiVals("EigenvaluesNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), *HEiValList);
#endif
//    for (int i = 0; i < SIZE; i++) {
//        delete hamilton1[i];
//    } delete[] hamilton1;
    }
}

/////////////////////////////// fixed magnetization blocks ///////////////////////////////

namespace magnetizationBlocks {
    void fillHamiltonBlock(const double &J1, const double &J2, const std::vector<int> &states, double **hamiltonBlock,
                      const int &N, const int &SIZE) { // remove continue; after else
        for (int i = 0; i < states.size(); i++) {
            int s = states.at(i);
            for (int n = 0; n < N / 2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    hamiltonBlock[i][i] += 0.25 * J1;
                } else {
                    hamiltonBlock[i][i] -= 0.25 * J1;
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int pos_d = findState(states, d);
                    hamiltonBlock[i][pos_d] = 0.5 * J1;
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    hamiltonBlock[i][i] += 0.25 * J2;
                } else {
                    hamiltonBlock[i][i] -= 0.25 * J2;
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    int pos_d = findState(states, d);
                    hamiltonBlock[i][pos_d] = 0.5 * J2;
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    hamiltonBlock[i][i] += 0.25 * J2;
                } else {
                    hamiltonBlock[i][i] -= 0.25 * J2;
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    int pos_d = findState(states, d);
                    hamiltonBlock[i][pos_d] = 0.5 * J2;
                }
            }
        }
    }

    void getEiValsFromBlock(const double &J1, const double &J2, const int &m,
                            std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks,
                            const int &N, const int &SIZE) {
        auto *states = new std::vector<int>;
        fillStates(states, m, N, SIZE);
        const int statesCount = (int) states->size();
        auto **hamiltonBlock = new double *[statesCount];
        for (int i = 0; i < statesCount; i++) {
            hamiltonBlock[i] = new double[statesCount];
            for (int j = 0; j < statesCount; j++) {
                hamiltonBlock[i][j] = 0.0;
            }
        }

        fillHamiltonBlock(J1, J2, *states, hamiltonBlock, N, SIZE);

        Eigen::MatrixXd H(statesCount, statesCount);
        for (int i = 0; i < statesCount; i++) {
            H.row(i) = Eigen::VectorXd::Map(&hamiltonBlock[i][0], statesCount);
        }

#if defined(showMatrix) || defined(saveMatrix)
        matrixBlocks->push_back(H);
#endif

        Eigen::EigenSolver<Eigen::MatrixXd> solver(H);
        const Eigen::VectorXcd &H1EiVal = solver.eigenvalues();
        for (std::complex<double> ev: H1EiVal) {
            HEiValList->push_back(ev);
        }

        states->clear();
        delete states;
        for (int i = 0; i < statesCount; i++) {
            delete hamiltonBlock[i];
        }
        delete[] hamiltonBlock;
    }

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE) {
        for (int m = 0; m <= N; m++) {
            getEiValsFromBlock(J1, J2, m, HEiValList, matrixBlocks, N, SIZE);
        }
    #if defined(showMatrix) || defined(saveMatrix)
        int offset_mag_blocks = 0;
        Eigen::MatrixXd H_mag_Block = Eigen::MatrixXd::Zero(SIZE, SIZE);
        for (const Eigen::MatrixXd &M: *matrixBlocks) {
            H_mag_Block.block(offset_mag_blocks, offset_mag_blocks, M.rows(), M.cols()) = M;
            offset_mag_blocks += (int) M.rows();
        }
    #endif
    #ifdef showMatrix
        std::cout << H_mag_Block << "\n";
    #endif
    #ifdef saveMatrix
        saveMatrixToFile(H_mag_Block, "HamiltonMagnetization.txt",
                         "Magnetisierungs Ansatz für N = " + std::to_string(N) +
                         "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2));
    #endif
    #if defined(showEigenvalues) || defined(saveEigenvalues)
        HEiValList->shrink_to_fit();
        // sort eigenvalues
        std::sort(HEiValList->begin(), HEiValList->end(),
                  [](const std::complex<double> &c1, const std::complex<double> &c2) {
                      return std::real(c1) < std::real(c2);
                  });
    #endif
    #ifdef showEigenvalues
        std::cout << "eigenvalues:\n";
        for (std::complex<double> ev: *HEiValList) {
            std::cout << ev << "\n";
        }
    #endif
    #ifdef saveEigenvalues
        saveComplexEiVals("EigenvaluesMagnetization.txt", "Magnetisierungs Ansatz für N = " + std::to_string(N)
                        + "\nJ1 = " + std::to_string(J1) + "\nJ2 = "
                        + std::to_string(J2), *HEiValList);
    #endif
    }
}

/////////////////////////////// momentum states ///////////////////////////////

namespace momentumStates {
    void fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                           const std::vector<int> &R_vals, std::complex<double> **hamiltonBlock, const int &N,
                           const int &SIZE) {
        for (int a = 0; a < states.size(); a++) {
            int s = states.at(a);
            for (int n = 0; n < N / 2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    hamiltonBlock[a][a] += std::complex<double>(0.25 * J1, 0.0);
                } else {
                    hamiltonBlock[a][a] -= std::complex<double>(0.25 * J1, 0.0);
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int r = 0, l = 0;
                    representative(d, &r, &l, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                        hamiltonBlock[a][b] +=
                                (std::complex<double>) 0.5 * J1 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) *
                                std::exp(numC);
                    }
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    hamiltonBlock[a][a] += std::complex<double>(0.25 * J2, 0.0);
                } else {
                    hamiltonBlock[a][a] -= std::complex<double>(0.25 * J2, 0.0);
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    int r = 0, l = 0;
                    representative(d, &r, &l, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                        hamiltonBlock[a][b] +=
                                (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) *
                                std::exp(numC);
                    }
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    hamiltonBlock[a][a] += std::complex<double>(0.25 * J2, 0.0);
                } else {
                    hamiltonBlock[a][a] -= std::complex<double>(0.25 * J2, 0.0);
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    int r = 0, l = 0;
                    representative(d, &r, &l, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                        hamiltonBlock[a][b] +=
                                (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) *
                                std::exp(numC);
                    }
                }
            }
        }
    }

    void momentumBlockSolver(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                             const std::vector<int> &R_vals, std::vector<std::complex<double>> *HEiValList,
                             std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE) {
        const int statesCount = (int) states.size();
        if (statesCount == 0) {
            return;
        }
        auto **hamiltonBlock = new std::complex<double> *[statesCount];
        for (int i = 0; i < statesCount; i++) {
            hamiltonBlock[i] = new std::complex<double>[statesCount];
            for (int j = 0; j < statesCount; j++) {
                hamiltonBlock[i][j] = 0.0;
            }
        }

        fillHamiltonBlock(J1, J2, k, states, R_vals, hamiltonBlock, N, SIZE);

        Eigen::MatrixXcd H(statesCount, statesCount);
        for (int i = 0; i < statesCount; i++) {
            H.row(i) = Eigen::VectorXcd::Map(&hamiltonBlock[i][0], statesCount);
        }

#if defined(showMatrix) || defined(saveMatrix)
        matrixBlocks->push_back(H);
#endif

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(H);
        const Eigen::VectorXcd &H1EiVal = solver.eigenvalues();
        for (std::complex<double> ev: H1EiVal) {
            HEiValList->push_back(ev);
        }

        for (int i = 0; i < statesCount; i++) {
            delete hamiltonBlock[i];
        }
        delete[] hamiltonBlock;
    }

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE) {

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<std::vector<std::vector<int>>> vec1(N + 1, std::vector<std::vector<int>>(N));
        std::vector<std::vector<std::vector<int>>> vec2(N + 1, std::vector<std::vector<int>>(N));
        auto *states = &vec1;
        auto *R_vals = &vec2;

        for (int s = 0; s < SIZE; s++) {
            int m = bitSum(s, N);
            for (int k = k_lower; k <= k_upper; k++) {
                int R = checkState(s, k, N);
                if (R >= 0) {
                    states->at(m).at(k - k_lower).push_back(s);
                    R_vals->at(m).at(k - k_lower).push_back(R);
                }
            }
        }

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                momentumBlockSolver(J1, J2, k, states->at(m).at(k - k_lower), R_vals->at(m).at(k - k_lower), HEiValList,
                                    matrixBlocks, N, SIZE);
            }
        }

    #if defined(showMatrix) || defined(saveMatrix)
        int offset_blocks = 0;
        Eigen::MatrixXcd H_moment_Block = Eigen::MatrixXcd::Zero(SIZE, SIZE);
        for (const Eigen::MatrixXcd &M: *matrixBlocks) {
            H_moment_Block.block(offset_blocks, offset_blocks, M.rows(), M.cols()) = M;
            offset_blocks += (int) M.rows();
        }
    #endif
    #ifdef showMatrix
        std::cout << H_moment_Block << "\n";
    #endif
    #ifdef saveMatrix
        saveComplexMatrixToFile(H_moment_Block, "HamiltonMomentumStates.txt",
                                "momentumStateAnsatz states Ansatz für N = " + std::to_string(N) +
                                "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2));
    #endif
    #if defined(showEigenvalues) || defined(saveEigenvalues)
        HEiValList->shrink_to_fit();
        // sort eigenvalues
        std::sort(HEiValList->begin(), HEiValList->end(),
                  [](const std::complex<double> &c1, const std::complex<double> &c2) {
                      return std::real(c1) < std::real(c2);
                  });
    #endif
    #ifdef showEigenvalues
        std::cout << "eigenvalues:\n";
        for (std::complex<double> ev: *HEiValList) {
            std::cout << ev << "\n";
        }
    #endif
    #ifdef saveEigenvalues
        saveComplexEiVals("EigenvaluesMomentumStates.txt",
                          "momentumStateAnsatz states Ansatz für N = " + std::to_string(N) +
                          "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *HEiValList);
    #endif
    }

}

/////////////////////////////// parity states (unfinished) ///////////////////////////////

#ifdef parityStateAnsatz
    void parityStateAnsatz(double J1, double J2, std::vector<std::complex<double>> *eiVals, std::vector<Eigen::MatrixXcd> *matrixBlocks) {
    int k_lower = -(N+2)/4+1;
    int k_upper = N/4;

    std::vector<std::vector<std::vector<int>>> vec(N+1, std::vector<std::vector<int>>(N));
    auto *states_mk = &vec;

    for (int s = 0; s < SIZE; s++) {
        int m = bitSum(s, N);
        for (int k = k_lower; k <= k_upper; k++) {
            int R = checkState(s, k, N);
            if (R >= 0) {
                states_mk->at(m).at(k-k_lower).push_back(s);
            }
        }
    }

    auto *states = new std::vector<int>;
    auto *R_vals = new std::vector<int>;
    auto *m_vals = new std::vector<int>;

    for (int m = 0; m <= N; m++) {
        for (int k = k_lower; k <= k_upper; k++) {
            for (int s : states_mk->at(m).at(k-k_lower)) {
                int R, m_cs;
                checkState(s, &R, &m_cs, k, N);
                for (int p : {-1, 1}) {
                    // p = +1 for k != 0, pi
                    if (k != 0 && k != k_upper && p == -1 ) {
                        continue;
                    }
                    for (int sigma : {-1, 1}) {
                        // sigma = +1 if k == 0 or k == k_upper (N/4)
                        if ((k == 0 || k == k_upper) && sigma == -1) {
                            continue;
                        }
                        if (m_cs != -1) {
                            std::complex<double> val  = (double) sigma * (double) p * std::cos(std::complex<double>(0, 4 * PI * (double) k * (double) m_cs / (double) N));
                            if (abs(std::complex<double>(1,0) + val) < 0.0001) {R = -1;}
                            if (sigma == -1 && abs(std::complex<double>(1,0) - val) > 0.0001) {R = -1;}
                        } if (R > 0) {
                            states->push_back(s);
                            R_vals->push_back(sigma * R);
                            m_vals->push_back(m_cs);
                        }

                    }
                    // ??????????
                }
            }


            parityBlockSolver(J1, J2, k, states->at(m).at(k-k_lower), R_vals->at(m).at(k-k_lower), eiVals, matrixBlocks);
        }
    }

    // sort eigenvalues
    std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
        return std::real(c1) < std::real(c2);
    });

#if defined(showMatrix) || defined(saveMatrix)
    int offset_blocks = 0;
    Eigen::MatrixXcd H_parity_Block = Eigen::MatrixXcd::Zero(SIZE, SIZE);
    for (const Eigen::MatrixXcd& M : *matrixBlocks) {
        H_parity_Block.block(offset_blocks, offset_blocks, M.rows(), M.cols()) = M;
        offset_blocks += M.rows();
    }
#endif
#ifdef showMatrix
    coutMutex.lock();
    std::cout << H_parity_Block << "\n";
    coutMutex.unlock();
#endif
#ifdef saveMatrix
    saveComplexMatrixToFile(H_parity_Block, "HamiltonParityStates.txt", "parityStateAnsatz states Ansatz für N = " + std::to_string(N) +
                                                               "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2));
#endif
#ifdef showEigenvalues
    coutMutex.lock();
    std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : *eiVals) {
        std::cout << ev << "\n";
    }
    coutMutex.unlock();
#endif
#ifdef saveEigenvalues
    saveComplexEiVals("EigenvaluesMParityStates.txt", "parityStateAnsatz states Ansatz für N = " + std::to_string(N) +
                                                       "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *eiVals);
#endif
}
#endif
