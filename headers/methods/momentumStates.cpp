#include "momentumStates.h"

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

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);
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

        std::vector<std::vector<std::vector<int>>> vec1(N + 1, std::vector<std::vector<int>>(N/2));
        std::vector<std::vector<std::vector<int>>> vec2(N + 1, std::vector<std::vector<int>>(N/2));
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
                                "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), N);
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
                          "momentum states Ansatz für N = " + std::to_string(N) +
                          "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *HEiValList, N);
#endif
    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "momentum states:..." << std::endl;

        auto *momentEiVals = new std::vector<std::complex<double>>;
        auto *matrixMomentBlocks = new std::vector<Eigen::MatrixXcd>;

        getEiVals(J1, J2, momentEiVals, matrixMomentBlocks, N, SIZE);

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        delete momentEiVals;
        delete matrixMomentBlocks;
    }

    void startSpecificHeat(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                           const double &END, const int &COUNT, const int &cores) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "specific heat (momentum states): calculating..." << std::endl;

        auto *momentEiVals = new std::vector<std::complex<double>>;
        auto *matrixMomentBlocks = new std::vector<Eigen::MatrixXcd>;

        getEiVals(J1, J2, momentEiVals, matrixMomentBlocks, N, SIZE);


        ///// specific /////

        auto *specificHeat_momentum = new std::vector<std::tuple<double, double>>;

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            specificHeat_momentum->push_back({current, getSpecificHeat(current, *momentEiVals, N)});
        }

        ///// save /////

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        std::string filenameSpecificHeat_C = "momentum_specific_heat.txt";
        std::string headerSpecificHeat_C = "N: " + std::to_string(N) + "\n"
                                           + "T START: " + std::to_string(START) + "\n"
                                           + "T END: " + std::to_string(END) + "\n"
                                           + "data-points: " + std::to_string(COUNT) + "\n"
                                           + "calculation time with " + std::to_string(cores) + " threads: " + formatTime(elapsed_seconds);

        std::string headerWithJSpecificHeat_C = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSpecificHeat_C;
        saveOutData(filenameSpecificHeat_C, headerWithJSpecificHeat_C, "J1/J2", "specific heat in J2", *specificHeat_momentum, N);

        std::cout << "\n";

        delete momentEiVals;
        delete matrixMomentBlocks;

    }

    void startDispersionPlot(const double &J1, const double &J2, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "dispersion (momentum states): calculating..." << std::endl;

        auto *momentEiVals = new std::vector<std::vector<std::complex<double>>>(N/2);
        auto *momentData = new std::vector<std::tuple<int, double>>;
        auto *matrixMomentBlocks = new std::vector<Eigen::MatrixXcd>;

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<std::vector<std::vector<int>>> vec1(N + 1, std::vector<std::vector<int>>(N/2));
        std::vector<std::vector<std::vector<int>>> vec2(N + 1, std::vector<std::vector<int>>(N/2));
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
                momentumBlockSolver(J1, J2, k, states->at(m).at(k-k_lower), R_vals->at(m).at(k-k_lower), &momentEiVals->at(k-k_lower),
                                    matrixMomentBlocks, N, SIZE);
            }
        }

        for (int k = k_lower; k <= k_upper; k++) {
            for (std::complex<double> ev : momentEiVals->at(k-k_lower)) {
                momentData->push_back({k, std::real(ev)});
            }
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        ///// save /////

        std::string filename = "momentum_energy_dispersion_J_const.txt";
        std::string header = "N: " + std::to_string(N);
        std::string headerWithJ = "J1/J2 = " + std::to_string(J1/J2) +"\n" + header;
        saveOutData(filename, headerWithJ, "k", "E in J2", *momentData, N);

        std::cout << "\n";

        delete momentEiVals;
        delete matrixMomentBlocks;

    }
}
