#include "momentumStates.h"

/////////////////////////////// momentum states ///////////////////////////////

namespace ED::momentumStates {

    void fillHamiltonBlock(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                           const std::vector<int> &R_vals, Eigen::MatrixXcd &hamiltonBlock, const int &N,
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
                    hamiltonBlock(a,a) += std::complex<double>(0.25 * J1, 0.0);
                } else {
                    hamiltonBlock(a,a) -= std::complex<double>(0.25 * J1, 0.0);
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int r = 0, l = 0;
                    representative(d, &r, &l, N);
                    int b = findState(states, r);
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
                    representative(d, &r, &l, N);
                    int b = findState(states, r);
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
                    representative(d, &r, &l, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                        hamiltonBlock(a,b) += (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a)
                                            / (double) R_vals.at(b)) * std::exp(numC);
                    }
                }
            }
        }
    }

    void momentumBlockSolver(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                             const std::vector<int> &R_vals, std::vector<std::complex<double>> &HEiValList,
                             std::vector<Eigen::MatrixXcd> &matrixBlocks, const int &N, const int &SIZE) {

        const int statesCount = (int) states.size();
        if (statesCount == 0) {
            return;
        }
        Eigen::MatrixXcd hamiltonBlock = Eigen::MatrixXcd::Zero(statesCount,statesCount);

//        std::cout << "M = " << bitSum(states.at(0), N) << ", k = " << k << std::endl;

        fillHamiltonBlock(J1, J2, k, states, R_vals, hamiltonBlock, N, SIZE);

        #if defined(showMatrix) || defined(saveMatrix)
            matrixBlocks.push_back(hamiltonBlock);
        #endif

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(hamiltonBlock);
        const Eigen::VectorXcd &H1EiVal = solver.eigenvalues();
        for (std::complex<double> ev: H1EiVal) {
            HEiValList.push_back(ev);
        }

    }

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> &HEiValList,
                   std::vector<Eigen::MatrixXcd> &matrixBlocks, const int &N, const int &SIZE) {

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<std::vector<std::vector<int>>> states(N + 1, std::vector<std::vector<int>>(N/2));
        std::vector<std::vector<std::vector<int>>> R_vals(N + 1, std::vector<std::vector<int>>(N/2));

        for (int s = 0; s < SIZE; s++) {
            int m = bitSum(s, N);
            for (int k = k_lower; k <= k_upper; k++) {
                int R = checkState(s, k, N);
                if (R >= 0) {
                    states.at(m).at(k - k_lower).push_back(s);
                    R_vals.at(m).at(k - k_lower).push_back(R);
                }
            }
        }

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                momentumBlockSolver(J1, J2, k, states.at(m).at(k - k_lower), R_vals.at(m).at(k - k_lower), HEiValList,
                                    matrixBlocks, N, SIZE);
            }
        }

        #if defined(showMatrix) || defined(saveMatrix)
            int offset_blocks = 0;
            Eigen::MatrixXcd H_moment_Block = Eigen::MatrixXcd::Zero(SIZE, SIZE);
            for (const Eigen::MatrixXcd &M: matrixBlocks) {
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
            HEiValList.shrink_to_fit();
            // sort eigenvalues
            std::sort(HEiValList.begin(), HEiValList.end(),
                      [](const std::complex<double> &c1, const std::complex<double> &c2) {
                          return std::real(c1) < std::real(c2);
                      });
        #endif
        #ifdef showEigenvalues
            std::cout << "eigenvalues:\n";
            for (std::complex<double> ev: HEiValList) {
                std::cout << ev << "\n";
            }
        #endif
        #ifdef saveEigenvalues
            saveComplexEiVals("EigenvaluesMomentumStates.txt",
                              "momentum states Ansatz für N = " + std::to_string(N) +
                              "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), HEiValList, N);
        #endif
    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "momentum states:..." << std::endl;

        std::vector<std::complex<double>> momentEiVals;
        std::vector<Eigen::MatrixXcd> matrixMomentBlocks;

        getEiVals(J1, J2, momentEiVals, matrixMomentBlocks, N, SIZE);

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

    }

    void startSpecificHeat(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                           const double &END, const int &COUNT) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "specific heat (momentum states): calculating..." << std::endl;

        std::vector<std::complex<double>> momentEiVals;
        std::vector<Eigen::MatrixXcd> matrixMomentBlocks;

        getEiVals(J1, J2, momentEiVals, matrixMomentBlocks, N, SIZE);


        ///// specific /////

        std::vector<std::tuple<double, double>> specificHeat_momentum;

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            specificHeat_momentum.emplace_back(current, getSpecificHeat(current, momentEiVals, N));
        }

        ///// save /////

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        std::string filenameSpecificHeat_C = "data_specific_heat_J_const.txt"; // momentum_specific_heat / data_specific_heat_J_const
        std::string headerSpecificHeat_C = "N: " + std::to_string(N) + "\n"
                                           + "T START: " + std::to_string(START) + "\n"
                                           + "T END: " + std::to_string(END) + "\n"
                                           + "data-points: " + std::to_string(COUNT) + "\n"
                                           + "calculation time: " + formatTime(elapsed_seconds);

        std::string headerWithJSpecificHeat_C = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSpecificHeat_C;
        saveOutData(filenameSpecificHeat_C, headerWithJSpecificHeat_C, "J1/J2", "specific heat in J2", specificHeat_momentum, N);
//        std::cout << "\n";

    }

    void startDispersionPlot(const double &J1, const double &J2, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "dispersion (momentum states): calculating..." << std::endl;

        std::vector<std::vector<std::complex<double>>> momentEiVals(N/2);
        std::vector<std::tuple<int, double>> momentData;
        std::vector<Eigen::MatrixXcd> matrixMomentBlocks;

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<std::vector<std::vector<int>>> states(N + 1, std::vector<std::vector<int>>(N/2));
        std::vector<std::vector<std::vector<int>>> R_vals(N + 1, std::vector<std::vector<int>>(N/2));

        for (int s = 0; s < SIZE; s++) {
            int m = bitSum(s, N);
            for (int k = k_lower; k <= k_upper; k++) {
                int R = checkState(s, k, N);
                if (R >= 0) {
                    states.at(m).at(k - k_lower).push_back(s);
                    R_vals.at(m).at(k - k_lower).push_back(R);
                }
            }
        }

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                momentumBlockSolver(J1, J2, k, states.at(m).at(k-k_lower), R_vals.at(m).at(k-k_lower), momentEiVals.at(k-k_lower),
                                    matrixMomentBlocks, N, SIZE);
            }
        }

        for (int k = k_lower; k <= k_upper; k++) {
            for (std::complex<double> ev : momentEiVals.at(k-k_lower)) {
                momentData.emplace_back(k, std::real(ev));
            }
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        ///// save /////

        std::string filename = "momentum_energy_dispersion_J_const.txt";
        std::string header = "N: " + std::to_string(N);
        std::string headerWithJ = "J1/J2 = " + std::to_string(J1/J2) +"\n" + header;
        saveOutData(filename, headerWithJ, "k", "E in J2", momentData, N);
//        std::cout << "\n";

    }

    void momentumBlockSolver_with_S(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                                        const std::vector<int> &R_vals, std::vector<std::tuple<std::complex<double>, int>> &data,
                                        const int &N, const int &SIZE) {

        const int statesCount = (int) states.size();
        if (statesCount == 0) {
            return;
        }
        Eigen::MatrixXcd hamiltonBlock = Eigen::MatrixXcd::Zero(statesCount,statesCount);

        fillHamiltonBlock(J1, J2, k, states, R_vals, hamiltonBlock, N, SIZE);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(hamiltonBlock);
        const Eigen::VectorXcd &H1EiVal = solver.eigenvalues();
        std::vector<std::complex<double>> HEiValList;
        for (std::complex<double> ev: H1EiVal) {
            HEiValList.push_back(ev);
        }

//        std::cout << "getting S2\n";
        Eigen::MatrixXcd S2 = spinMatrixMomentum(N, k, states, R_vals);
//        std::cout << "getting U\n";
        const Eigen::MatrixXcd& U = solver.eigenvectors();
//        std::cout << "getting U_inv_S2_U\n";
        Eigen::MatrixXcd U_inv_S2_U = U.adjoint() * S2 * U;

//        std::cout << "filling data in momentumBlockSolver_with_S\n";
//        std::cout << "HEiValList size: " << HEiValList.size();
        for (int i = 0; i < HEiValList.size(); i++) {
            double S_elem = std::real(U_inv_S2_U(i, i));
            double S = - 0.5 + std::sqrt(0.25 + S_elem);
            if (S < EPSILON) { S = 0.0;}
            data.emplace_back(HEiValList.at(i), S);
//            std::cout << S << std::endl;
        }
//        std::cout << ", data size: " << data.size() << std::endl;


    }

    void getEiValsZeroBlock(const double &J1, const double &J2, std::vector<std::vector<std::complex<double>>> &eiVals,
                            std::vector<Eigen::MatrixXcd> &matrixBlockU, std::vector<Eigen::MatrixXcd> &matrixBlockS2,
                            const int &N, const int &SIZE) {

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<int> states_m;
        std::vector<std::vector<int>> states(N/2);
        std::vector<std::vector<int>> R_vals(N/2);
        fillStates(&states_m, N/2, N, SIZE);

        for (int s : states_m) {
            for (int k = k_lower; k <= k_upper; k++) {
                int R = checkState(s, k, N);
                if (R >= 0) {
                    states.at(k - k_lower).push_back(s);
                    R_vals.at(k - k_lower).push_back(R);
                }
            }
        }

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                momentumBlockSolver_withMatrix(J1, J2, k, states.at(k - k_lower), R_vals.at(k - k_lower), eiVals,
                                    matrixBlockU, matrixBlockS2, N, SIZE);
            }
        }

    }

    void momentumBlockSolver_withMatrix(const double &J1, const double &J2, const int &k, const std::vector<int> &states,
                                        const std::vector<int> &R_vals, std::vector<std::vector<std::complex<double>>> &eiVals,
                                        std::vector<Eigen::MatrixXcd> &matrixBlockU, std::vector<Eigen::MatrixXcd> &matrixBlockS2,
                                        const int &N, const int &SIZE) {

        const int statesCount = (int) states.size();
        if (statesCount == 0) {
            return;
        }
        Eigen::MatrixXcd hamiltonBlock = Eigen::MatrixXcd::Zero(statesCount,statesCount);

        fillHamiltonBlock(J1, J2, k, states, R_vals, hamiltonBlock, N, SIZE);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(hamiltonBlock);
        matrixBlockU.push_back(solver.eigenvectors());
        const Eigen::VectorXcd &H1EiVal = solver.eigenvalues();
        std::vector<std::complex<double>> HEiValList;
        for (std::complex<double> ev: H1EiVal) {
            HEiValList.push_back(ev);
        }
        HEiValList.shrink_to_fit();
        eiVals.push_back(HEiValList);

        Eigen::MatrixXcd S2 = spinMatrixMomentum(N, k, states, R_vals);
        matrixBlockS2.push_back(S2);

    }

    void getEiValsMagBlock_with_index(const double &J1, const double &J2, std::vector<std::tuple<std::complex<double>, int, int>> &data,
                            const int &N, const int &SIZE, const int &mag) {

        std::vector<Eigen::MatrixXcd> matrixBlockU;

        int k_lower = 0;//-(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<int> states_m;
        std::vector<std::vector<int>> states(N/2);
        std::vector<std::vector<int>> R_vals(N/2);
        fillStates(&states_m, mag, N, SIZE);

        for (int s : states_m) {
            for (int k = k_lower; k <= k_upper; k++) {
                int R = checkState(s, k, N);
                if (R >= 0) {
                    states.at(k - k_lower).push_back(s);
                    R_vals.at(k - k_lower).push_back(R);
                }
            }
        }

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                std::vector<std::tuple<std::complex<double>, int>> eiVals;
//                std::cout << "calling momentumBlockSolver_with_S\n";
                momentumBlockSolver_with_S(J1, J2, k, states.at(k - k_lower), R_vals.at(k - k_lower), eiVals, N, SIZE);
//                std::cout << "eivals size: " << eiVals.size();
//                std::cout << "filling data in getEiValsMagBlock_with_index\n";
                for (auto & eiVal : eiVals) {
                    data.emplace_back(std::get<0>(eiVal), k, std::get<1>(eiVal));
                }
//                std::cout << ", data size: " << data.size() << std::endl;
//                std::cout << "filled data in getEiValsMagBlock_with_index\n";
            }
        }

    }

    void getEiValsMagBlock(const double &J1, const double &J2, std::vector<std::complex<double>> &eiVals,
                                  const int &N, const int &SIZE, const int &mag) {

        std::vector<Eigen::MatrixXcd> matrixBlockU;

        int k_lower = -(N + 2) / 4 + 1;
        int k_upper = N / 4;

        std::vector<int> states_m;
        std::vector<std::vector<int>> states(N/2);
        std::vector<std::vector<int>> R_vals(N/2);
        fillStates(&states_m, mag, N, SIZE);

        for (int s : states_m) {
            for (int k = k_lower; k <= k_upper; k++) {
                int R = checkState(s, k, N);
                if (R >= 0) {
                    states.at(k - k_lower).push_back(s);
                    R_vals.at(k - k_lower).push_back(R);
                }
            }
        }

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                momentumBlockSolver(J1, J2, k, states.at(k - k_lower), R_vals.at(k - k_lower), eiVals,
                                    matrixBlockU, N, SIZE);
            }
        }

    }

    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "susceptibility (momentum states): calculating..." << std::endl;

        std::vector<std::vector<std::complex<double>>> eiVals;
        std::vector<Eigen::MatrixXcd> matrixBlockU;
        std::vector<Eigen::MatrixXcd> matrixBlockS2;
        getEiValsZeroBlock(J1, J2, eiVals, matrixBlockU, matrixBlockS2, N, SIZE);

        ///// susceptibility /////

        std::vector<std::tuple<double, double>> susceptibility_magnetization;

        std::vector<Eigen::MatrixXcd> Blocks_U_inv_S2_U;
        //std::cout << "M:\n";
        for(int i = 0; i < matrixBlockU.size(); i++) {
            Eigen::MatrixXcd M = matrixBlockU.at(i).adjoint() * matrixBlockS2.at(i) * matrixBlockU.at(i);
            Blocks_U_inv_S2_U.push_back(M);
//            for (int j = 0; j < M.rows(); j++) {
//                std::cout << M(j,j) << std::endl;
//            }
            //std::cout << M << std::endl;
        }

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            //current_beta = 1 / current_beta;
            susceptibility_magnetization.emplace_back(current, getSusceptibilityDegeneracy(current, Blocks_U_inv_S2_U, eiVals, N));
        }

//        for (std::tuple<double, double> t : susceptibility_magnetization) {
//            std::cout << std::get<1>(t) << "\n";
//        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        ///// save /////

        std::string filenameSusceptibility_X = "momentum_susceptibility_J_const.txt"; // momentum_susceptibility_J_const.txt / data_susceptibility_J_const
        std::string headerSusceptibility_X = "N: " + std::to_string(N) + "\n"
                                             + "T START: " + std::to_string(START) + "\n"
                                             + "T END: " + std::to_string(END) + "\n"
                                             + "data-points: " + std::to_string(COUNT) + "\n"
                                             + "calculation time: " + formatTime(elapsed_seconds);

        std::string headerWithJSusceptibility_X = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSusceptibility_X;
        saveOutData(filenameSusceptibility_X, headerWithJSusceptibility_X, "J1/J2", "specific heat in J2", susceptibility_magnetization, N);

        std::cout << "\n";

    }

}
