#include "magnetizationBlocks.h"

/////////////////////////////// fixed magnetization blocks ///////////////////////////////

namespace magnetizationBlocks {
    void fillHamiltonBlock(const double &J1, const double &J2, const std::vector<int> &states, double **hamiltonBlock,
                           const int &N) {

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

        fillHamiltonBlock(J1, J2, *states, hamiltonBlock, N);

        Eigen::MatrixXd H(statesCount, statesCount);
        for (int i = 0; i < statesCount; i++) {
            H.row(i) = Eigen::VectorXd::Map(&hamiltonBlock[i][0], statesCount);
        }

#if defined(showMatrix) || defined(saveMatrix)
        matrixBlocks->push_back(H);
#endif

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
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
        saveComplexEiVals("EigenvaluesMagnetization.txt", "Magnetisierungs Ansatz für N = " + std::to_string(N)
                        + "\nJ1 = " + std::to_string(J1) + "\nJ2 = "
                        + std::to_string(J2), *HEiValList, N);
#endif
    }

    void getEiValsZeroBlock(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                            Eigen::MatrixXcd &matrixBlockU, const std::vector<int> &states, const int &N) {

        const int statesCount = (int) states.size();
        auto **hamiltonBlock = new double *[statesCount];
        for (int i = 0; i < statesCount; i++) {
            hamiltonBlock[i] = new double[statesCount];
            for (int j = 0; j < statesCount; j++) {
                hamiltonBlock[i][j] = 0.0;
            }
        }

        fillHamiltonBlock(J1, J2, states, hamiltonBlock, N);

        Eigen::MatrixXd H(statesCount, statesCount);
        for (int i = 0; i < statesCount; i++) {
            H.row(i) = Eigen::VectorXd::Map(&hamiltonBlock[i][0], statesCount);
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
        const Eigen::VectorXcd &H1EiVal = solver.eigenvalues();
        for (std::complex<double> ev: H1EiVal) {
            HEiValList->push_back(ev);
        }

        matrixBlockU = solver.eigenvectors();

        for (int i = 0; i < statesCount; i++) {
            delete hamiltonBlock[i];
        }
        delete[] hamiltonBlock;

        HEiValList->shrink_to_fit();

    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "block diagonale m_z:..." << std::endl;
        auto *EiVals_m = new std::vector<std::complex<double>>;
        auto *matrixBlocks_m = new std::vector<Eigen::MatrixXd>;
        magnetizationBlocks::getEiVals(J1, J2, EiVals_m, matrixBlocks_m, N, SIZE);

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        delete matrixBlocks_m;
    }

    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "susceptibility (m_z blocks): calculating..." << std::endl;

        auto *states = new std::vector<int>;
        fillStates(states, N/2, N, SIZE);
        const int statesCount = (int) states->size();

        auto *eiVals = new std::vector<std::complex<double>>;
        Eigen::MatrixXcd matrixBlockU;
        magnetizationBlocks::getEiValsZeroBlock(J1, J2, eiVals, matrixBlockU, *states, N);

        for (std::complex<double> ev : *eiVals) {
            std::cout << ev << "\n";
        }

        ///// susceptibility /////

        auto *susceptibility_magnetization = new std::vector<std::tuple<double, double>>;

        Eigen::MatrixXd S2 = spinMatrix(N, *states);
        Eigen::MatrixXcd Matrix_U_inv_S2_U = Eigen::MatrixXcd::Zero(statesCount, statesCount);
        Matrix_U_inv_S2_U = matrixBlockU.adjoint() * S2 * matrixBlockU;

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            //current_beta = 1 / current_beta;
            susceptibility_magnetization->push_back({current, getSusceptibilityDegeneracy(current, Matrix_U_inv_S2_U, *eiVals, N)});
        }


        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        ///// save /////

        std::string filenameSusceptibility_X = "magnetization_susceptibility_J_const.txt";
        std::string headerSusceptibility_X = "N: " + std::to_string(N) + "\n"
                                             + "T START: " + std::to_string(START) + "\n"
                                             + "T END: " + std::to_string(END) + "\n"
                                             + "data-points: " + std::to_string(COUNT) + "\n"
                                             + "calculation time: " + formatTime(elapsed_seconds);

        std::string headerWithJSusceptibility_X = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSusceptibility_X;
        saveOutData(filenameSusceptibility_X, headerWithJSusceptibility_X, "J1/J2", "specific heat in J2", *susceptibility_magnetization, N);

        std::cout << "\n";

    }
}
