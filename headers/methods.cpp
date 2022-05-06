#include "methods.h"
#include "helpers.h"

/////////////////////////////// multi-threading ///////////////////////////////

namespace multi {
    void get_DeltaE_CT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                             double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                             const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        while (true) {

            int p = (int) ( (float) pos / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) pos / (float) COUNT * 100.0 ) << "% J1/J2 = " << J << " (" << pos << "/" << COUNT << ")     ";
            coutMutex.unlock();

            auto *eiVals = new std::vector<std::complex<double>>;
            auto *matrixBlocks = new std::vector<Eigen::MatrixXcd>;
            //auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;

            //naiv::getEiVals(J, 1.0, eiVals, N, SIZE);
            //magnetizationBlocks::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);
            momentumStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);

            // sort eigenvalues
            eiVals->shrink_to_fit();
            std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
                return std::real(c1) < std::real(c2);
            });

            // Delta E
            double E0 = std::real(eiVals->at(0));
            double E1 = std::real(eiVals->at(1));

            // write data
            nextJMutex.lock();
            outDataDeltaE->push_back({J, E1 - E0});
            outDataSpecificHeat_C->push_back({J, getSpecificHeat(T, *eiVals, N)});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                delete eiVals;
                delete matrixBlocks;
                break;
            } else {
                J = START + (END - START) * pos / COUNT;
            }
            // clean up
            eiVals->clear();
            matrixBlocks->clear();

        }

    }

    void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END, const unsigned int *cpu_cnt,
                               int *cores, const double &T, const int &N, const int &SIZE) {

        const clock_t begin_time = clock();

        std::cout << "\ncalculating:..." << std::endl;

        auto *outDataDeltaE = new std::vector<std::tuple<double, double>>;
        auto *outDataSpecificHeat_C = new std::vector<std::tuple<double, double>>;

        if (COUNT < *cpu_cnt) {
            *cores = COUNT;
        }
        std::thread Threads[*cores];

        CURRENT += *cores;

        for (int i = 0; i < *cores; i++) {
            Threads[i] = std::thread(get_DeltaE_CT_const, START + (END - START) * i / COUNT, i + 1, outDataDeltaE, T,
                                     outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
        }

        for (int i = 0; i < *cores; i++) {
            Threads[i].join();
        }

        auto time = float(clock () - begin_time) /  CLOCKS_PER_SEC;
        std::cout << "\n" << "calculations done; this took: " << time << " seconds\n\n";

        // sort data-points
        std::sort(outDataDeltaE->begin(), outDataDeltaE->end(), [](
                const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
                    return std::get<0>(a) < std::get<0>(b);
        });
        std::sort(outDataSpecificHeat_C->begin(), outDataSpecificHeat_C->end(), [](
                const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
                    return std::get<0>(a) < std::get<0>(b);
        });

        std::string filenameDeltaE = "data_delta_E.txt";
        std::string filenameSpecificHeat_C = "data_specific_heat.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(*cores) + " threads: " + std::to_string(time) + " seconds";

        std::string headerWithBeta = "T = " + std::to_string(T) + "\n" + header;

        saveOutData(filenameDeltaE, header, "J1/J2", "Delta E in J2", *outDataDeltaE, N);
        saveOutData(filenameSpecificHeat_C, headerWithBeta, "J1/J2", "specific heat in J2", *outDataSpecificHeat_C, N);

        delete outDataDeltaE;
        delete outDataSpecificHeat_C;

    }

    void get_XT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                      const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE, const double &T) {

        while (true) {

            int p = (int) ( (float) pos / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) pos / (float) COUNT * 100.0 ) << "% J1/J2 = " << J << " (" << pos << "/" << COUNT << ")     ";
            coutMutex.unlock();


            auto *states = new std::vector<int>;
            auto *eiVals = new std::vector<std::complex<double>>;
            Eigen::MatrixXcd matrixBlockU;
            magnetizationBlocks::getEiValsZeroBlock(J, 1.0, eiVals, matrixBlockU, states, N, SIZE);

            ///// susceptibility /////

            Eigen::MatrixXd S2 = spinMatrix(N, *states);
            Eigen::MatrixXcd Matrix_U_inv_S2_U = Eigen::MatrixXcd::Zero(SIZE, SIZE);
            Matrix_U_inv_S2_U = matrixBlockU.adjoint() * S2 * matrixBlockU;

            double susceptibility = getSusceptibilityDegeneracy(T, Matrix_U_inv_S2_U, *eiVals, N);

            // write data
            nextJMutex.lock();
            outDataMagneticSusceptibility_X->push_back({J, susceptibility});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                delete eiVals;
                delete states;
                break;
            } else {
                J = START + (END - START) * pos / COUNT;
            }
            // clean up
            eiVals->clear();
            states->clear();

        }

    }

    void start_XT_const(const int &COUNT, const double &START, const double &END, const unsigned int *cpu_cnt,
                        int *cores, const double &T, const int &N, const int &SIZE) {

        const clock_t begin_time = clock();

        std::cout << "\ncalculating:..." << std::endl;

        auto *outDataMagneticSusceptibility_X = new std::vector<std::tuple<double, double>>;

        if (COUNT < *cpu_cnt) {
            *cores = COUNT;
        }
        std::thread Threads[*cores];

        CURRENT = 1 + *cores;

        for (int i = 0; i < *cores; i++) {
            Threads[i] = std::thread(get_XT_const, START + (END - START) * i / COUNT, i + 1, outDataMagneticSusceptibility_X, COUNT, START, END, N, SIZE, T);
        }

        for (int i = 0; i < *cores; i++) {
            Threads[i].join();
        }

        auto time = float(clock () - begin_time) /  CLOCKS_PER_SEC;
        std::cout << "\n" << "calculations done; this took: " << time << " seconds\n\n";

        // sort data-points
        std::sort(outDataMagneticSusceptibility_X->begin(), outDataMagneticSusceptibility_X->end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::string filenameMagneticSusceptibility_X = "data_magnetic_susceptibility_T_const.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(*cores) + " threads: " + std::to_string(time) + " seconds";

        std::string headerWithBeta = "T = " + std::to_string(T) + "\n" + header;

        saveOutData(filenameMagneticSusceptibility_X, headerWithBeta, "J1/J2", "magnetic susceptibility in J2", *outDataMagneticSusceptibility_X, N);

        delete outDataMagneticSusceptibility_X;
    }

}

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
                   const int &N, const int &SIZE, Eigen::MatrixXcd &Matrix_U) {
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
        saveHamilton(hamilton1, "HamiltonNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), SIZE, N);
#endif
        //std::cout << "solving...\n";
        Eigen::EigenSolver<Eigen::MatrixXd> solver(H1);
        const Eigen::VectorXcd &H1EiVal = solver.eigenvalues();
        for (std::complex<double> ev : H1EiVal) {
            HEiValList->push_back(ev);
        }

        Matrix_U = solver.eigenvectors();

        HEiValList->shrink_to_fit();
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
        saveComplexEiVals("EigenvaluesNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), *HEiValList, N);
#endif
//    for (int i = 0; i < SIZE; i++) {
//        delete hamilton1[i];
//    } delete[] hamilton1;
    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
               const double &END, const int &COUNT, const int &cores) {
        const clock_t begin_time_NAIV = clock();

        std::cout << "\nnaiver Ansatz:..." << std::endl;
        auto *H_naiv_EiVals = new std::vector<std::complex<double>>;
        Eigen::MatrixXcd Matrix_U(SIZE, SIZE);
        getEiVals(J1, J2, H_naiv_EiVals, N, SIZE, Matrix_U);

        ///// susceptibility /////

        auto *susceptibility_naiv = new std::vector<std::tuple<double, double>>;

        Eigen::MatrixXd S2 = spinMatrix(N, SIZE);
        Eigen::MatrixXcd Matrix_U_inv_S2_U = Eigen::MatrixXcd::Zero(SIZE, SIZE);
        Matrix_U_inv_S2_U = Matrix_U.adjoint() * S2 * Matrix_U;

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            //current_beta = 1 / current_beta;
            susceptibility_naiv->push_back({current,
                                            getSusceptibility(current, Matrix_U_inv_S2_U, *H_naiv_EiVals, N)});
        }

        auto time_NAIV = float(clock () - begin_time_NAIV) /  CLOCKS_PER_SEC;
        std::cout << "calculations done; this took: " << time_NAIV << " seconds\n";

        ///// save /////

        std::string filenameSusceptibility_X = "naiv_susceptibility.txt";
        std::string headerSusceptibility_X = "N: " + std::to_string(N) + "\n"
                                             + "T START: " + std::to_string(START) + "\n"
                                             + "T END: " + std::to_string(END) + "\n"
                                             + "data-points: " + std::to_string(COUNT) + "\n"
                                             + "calculation time with " + std::to_string(cores) + " threads: " + std::to_string(time_NAIV) + " seconds";

        std::string headerWithJSusceptibility_X = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSusceptibility_X;
        saveOutData(filenameSusceptibility_X, headerWithJSusceptibility_X, "J1/J2", "specific heat in J2", *susceptibility_naiv, N);
        std::cout << "\n";
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
                   Eigen::MatrixXcd &matrixBlockU, std::vector<int> *states, const int &N, const int &SIZE) {

        fillStates(states, N/2, N, SIZE);
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

        Eigen::EigenSolver<Eigen::MatrixXd> solver(H);
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
        // sort eigenvalues
//        std::sort(HEiValList->begin(), HEiValList->end(),
//                  [](const std::complex<double> &c1, const std::complex<double> &c2) {
//                      return std::real(c1) < std::real(c2);
//                  });
    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {
        const clock_t begin_time_MAGNETIZATION = clock();

        std::cout << "\nblock diagonale m_z:..." << std::endl;
        //std::vector<std::complex<double>> EiVals_m;
        auto *EiVals_m = new std::vector<std::complex<double>>;
        auto *matrixBlocks_m = new std::vector<Eigen::MatrixXd>;
        magnetizationBlocks::getEiVals(J1, J2, EiVals_m, matrixBlocks_m, N, SIZE);

        auto time_MAGNETIZATION = float(clock () - begin_time_MAGNETIZATION) /  CLOCKS_PER_SEC;
        std::cout << "calculations done; this took: " << time_MAGNETIZATION << " seconds\n";
        delete matrixBlocks_m;
    }

    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT, const int &cores) {

        const clock_t begin_time_MAGNETIZATION = clock();
        std::cout << "\nblock diagonale m_z:..." << std::endl;

        auto *states = new std::vector<int>;
        auto *eiVals = new std::vector<std::complex<double>>;
        Eigen::MatrixXcd matrixBlockU;
        magnetizationBlocks::getEiValsZeroBlock(J1, J2, eiVals, matrixBlockU, states, N, SIZE);

        ///// susceptibility /////

        auto *susceptibility_magnetization = new std::vector<std::tuple<double, double>>;

        Eigen::MatrixXd S2 = spinMatrix(N, *states);
        Eigen::MatrixXcd Matrix_U_inv_S2_U = Eigen::MatrixXcd::Zero(SIZE, SIZE);
        Matrix_U_inv_S2_U = matrixBlockU.adjoint() * S2 * matrixBlockU;

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            //current_beta = 1 / current_beta;
            susceptibility_magnetization->push_back({current,
                                                     getSusceptibilityDegeneracy(current, Matrix_U_inv_S2_U, *eiVals, N)});
        }


        auto time_MAGNETIZATION = float(clock () - begin_time_MAGNETIZATION) /  CLOCKS_PER_SEC;
        std::cout << "calculations done; this took: " << time_MAGNETIZATION << " seconds\n";

        ///// save /////

        std::string filenameSusceptibility_X = "magnetization_susceptibility_J_const.txt";
        std::string headerSusceptibility_X = "N: " + std::to_string(N) + "\n"
                                             + "T START: " + std::to_string(START) + "\n"
                                             + "T END: " + std::to_string(END) + "\n"
                                             + "data-points: " + std::to_string(COUNT) + "\n"
                                             + "calculation time with " + std::to_string(cores) + " threads: " + std::to_string(time_MAGNETIZATION) + " seconds";

        std::string headerWithJSusceptibility_X = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSusceptibility_X;
        saveOutData(filenameSusceptibility_X, headerWithJSusceptibility_X, "J1/J2", "specific heat in J2", *susceptibility_magnetization, N);
        std::cout << "\n";
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
                          "momentumStateAnsatz states Ansatz für N = " + std::to_string(N) +
                          "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *HEiValList, N);
    #endif
    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {

        const clock_t begin_time_MOMENTUM = clock();

        std::cout << "\nmomentum states:..." << std::endl;

        auto *momentEiVals = new std::vector<std::complex<double>>;
        auto *matrixMomentBlocks = new std::vector<Eigen::MatrixXcd>;

        getEiVals(J1, J2, momentEiVals, matrixMomentBlocks, N, SIZE);

        auto time_MOMENTUM = float(clock () - begin_time_MOMENTUM) /  CLOCKS_PER_SEC;
        std::cout << "calculations done; this took: " << time_MOMENTUM << " seconds\n";

        std::cout << "\n";

        delete momentEiVals;
        delete matrixMomentBlocks;
    }

    void startSpecificHeat(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                           const double &END, const int &COUNT, const int &cores) {

        const clock_t begin_time_MOMENTUM = clock();

        std::cout << "\nmomentum states:..." << std::endl;

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

        auto time_MOMENTUM = float(clock () - begin_time_MOMENTUM) /  CLOCKS_PER_SEC;
        std::cout << "calculations done; this took: " << time_MOMENTUM << " seconds\n";

        std::string filenameSpecificHeat_C = "momentum_specific_heat.txt";
        std::string headerSpecificHeat_C = "N: " + std::to_string(N) + "\n"
                                           + "T START: " + std::to_string(START) + "\n"
                                           + "T END: " + std::to_string(END) + "\n"
                                           + "data-points: " + std::to_string(COUNT) + "\n"
                                           + "calculation time with " + std::to_string(cores) + " threads: " + std::to_string(time_MOMENTUM) + " seconds";

        std::string headerWithJSpecificHeat_C = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSpecificHeat_C;
        saveOutData(filenameSpecificHeat_C, headerWithJSpecificHeat_C, "J1/J2", "specific heat in J2", *specificHeat_momentum, N);

        std::cout << "\n";

        delete momentEiVals;
        delete matrixMomentBlocks;
    }
}

/////////////////////////////// parity states (unfinished) ///////////////////////////////

namespace parityStates {
    void getEiVals(double J1, double J2, std::vector<std::complex<double>> *eiVals,
                           std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE) {

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


            //parityBlockSolver(J1, J2, k, states->at(m).at(k-k_lower), R_vals->at(m).at(k-k_lower), eiVals, matrixBlocks);
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
    std::cout << H_parity_Block << "\n";
#endif
#ifdef saveMatrix
    saveComplexMatrixToFile(H_parity_Block, "HamiltonParityStates.txt", "parityStateAnsatz states Ansatz für N = " + std::to_string(N) +
                                                               "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), N);
#endif
#ifdef showEigenvalues
    std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : *eiVals) {
        std::cout << ev << "\n";
    }
#endif
#ifdef saveEigenvalues
    saveComplexEiVals("EigenvaluesMParityStates.txt", "parityStateAnsatz states Ansatz für N = " + std::to_string(N) +
                                                       "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *eiVals, N);
#endif
}

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {
        const clock_t begin_time_PARITY = clock();

        std::cout << "\nparity states:..." << std::endl;

        auto *parityEiVals = new std::vector<std::complex<double>>;
        auto *matrixParityBlocks = new std::vector<Eigen::MatrixXcd>;

        getEiVals(J1, J2, parityEiVals, matrixParityBlocks, N, SIZE);

        auto time_PARITY = float(clock () - begin_time_PARITY) /  CLOCKS_PER_SEC;
        std::cout << "calculations done; this took: " << time_PARITY << " seconds\n";
        delete parityEiVals;
        delete matrixParityBlocks;
    }
}