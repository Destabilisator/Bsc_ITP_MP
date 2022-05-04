#include "main.h"

//// methods ////
//#define naiv
//#define magnetization
#define momentum
//#define parity

#define multiCalc

///// output, turn off during multithreading /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

/////////////////////////////// naiver Ansatz ///////////////////////////////

#ifdef naiv
void fillHamiltonNaiv(double** hamilton, double J1, double J2) {
    for (int s = 0; s <= SIZE -1; s++) {
        for (int n = 0; n < N/2; n++) {
            // declaring indices
            int j_0 = 2 * n;
            int j_1 = (j_0+1) % N;
            int j_2 = (j_0+2) % N;
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

void naiverAnatz(double J1, double J2, std::vector<std::complex<double>> *HEiValList) {
    static auto **hamilton1 = new double*[SIZE];
    for (int i = 0; i < SIZE; i++) {
        hamilton1[i] = new double[SIZE];
        for (int j = 0; j < SIZE; j++) {
            hamilton1[i][j] = 0.0;
        }
    }
    fillHamiltonNaiv(hamilton1, J1, J2);
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
#endif

/////////////////////////////// fixed magnetization states ///////////////////////////////

#ifdef magnetization
void fillHamiltonBlock(double J1, double J2, const std::vector<int>& states, double **hamiltonBlock) { // remove continue; after else
    for (int i = 0; i < states.size(); i++) {
        int s = states.at(i);
        for (int n = 0; n < N/2; n++) {
            // declaring indices
            int j_0 = 2 * n;
            int j_1 = (j_0+1) % N;
            int j_2 = (j_0+2) % N;
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

void magBlock_getEiVal(double J1, double J2, int m, std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks) {
    auto *states = new std::vector<int>;
    fillStates(states, m, N, SIZE);
    const int statesCount = states->size();
    auto **hamiltonBlock = new double*[statesCount];
    for (int i = 0; i < statesCount; i++) {
        hamiltonBlock[i] = new double[statesCount];
        for (int j = 0; j < statesCount; j++) {
            hamiltonBlock[i][j] = 0.0;
        }
    }

    fillHamiltonBlock(J1, J2, *states, hamiltonBlock);

    Eigen::MatrixXd H(statesCount, statesCount);
    for (int i = 0; i < statesCount; i++) {
        H.row(i) = Eigen::VectorXd::Map(&hamiltonBlock[i][0], statesCount);
    }

#if defined(showMatrix) || defined(saveMatrix)
    matrixBlocks->push_back(H);
#endif

    Eigen::EigenSolver<Eigen::MatrixXd> solver(H);
    const Eigen::VectorXcd& H1EiVal = solver.eigenvalues();
    for (std::complex<double> ev : H1EiVal) {
        HEiValList->push_back(ev);
    }

    states->clear();
    delete states;
    for (int i = 0; i < statesCount; i++) {
        delete hamiltonBlock[i];
    }
    delete[] hamiltonBlock;
}

void magnetisierungsAnsatz(double J1, double J2, std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks) {
    for (int m = 0; m <= N; m++) {
        magBlock_getEiVal(J1, J2, m, HEiValList, matrixBlocks);
    }
#if defined(showMatrix) || defined(saveMatrix)
    int offset_mag_blocks = 0;
    Eigen::MatrixXd H_mag_Block = Eigen::MatrixXd::Zero(SIZE, SIZE);
    for (const Eigen::MatrixXd& M : *matrixBlocks) {
        H_mag_Block.block(offset_mag_blocks, offset_mag_blocks, M.rows(), M.cols()) = M;
        offset_mag_blocks += M.rows();
    }
#endif
#ifdef showMatrix
    coutMutex.lock();
    std::cout << H_mag_Block << "\n";
    coutMutex.unlock();
#endif
#ifdef saveMatrix
    saveMatrixToFile(H_mag_Block, "HamiltonMagnetization.txt", "Magnetisierungs Ansatz für N = " + std::to_string(N) +
                                                         "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2));
#endif
#if defined(showEigenvalues) || defined(saveEigenvalues)
    HEiValList->shrink_to_fit();
    // sort eigenvalues
    std::sort(HEiValList->begin(), HEiValList->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
        return std::real(c1) < std::real(c2);
    });
#endif
#ifdef showEigenvalues
    coutMutex.lock();
    std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : *HEiValList) {
        std::cout << ev << "\n";
    }
    coutMutex.unlock();
#endif
#ifdef saveEigenvalues
    saveComplexEiVals("EigenvaluesMagnetization.txt", "Magnetisierungs Ansatz für N = " + std::to_string(N) +
                      "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *HEiValList);
#endif
}
#endif

/////////////////////////////// momentum states ///////////////////////////////

#if defined(momentum) || defined(multiCalc)
void fillHamiltonMomentumBlock(double J1, double J2, int k,const std::vector<int> &states, const std::vector<int> &R_vals, std::complex<double> **hamiltonBlock) {
    for (int a = 0; a < states.size(); a++) {
        int s = states.at(a);
        for (int n = 0; n < N/2; n++) {
            // declaring indices
            int j_0 = 2 * n;
            int j_1 = (j_0+1) % N;
            int j_2 = (j_0+2) % N;
            // applying H to state s
            if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                hamiltonBlock[a][a] += std::complex<double> (0.25 * J1, 0.0);
            } else {
                hamiltonBlock[a][a] -= std::complex<double> (0.25 * J1, 0.0);
                int d = s ^ (1 << j_0) ^ (1 << j_2);
                int r = 0, l = 0;
                representative(d, &r, &l, N);
                int b = findState(states, r);
                if (b >= 0) {
                    std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N );
                    hamiltonBlock[a][b] += (std::complex<double>) 0.5 * J1 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * std::exp(numC);
                }
            }
            if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                hamiltonBlock[a][a] += std::complex<double> (0.25 * J2, 0.0);
            } else {
                hamiltonBlock[a][a] -= std::complex<double> (0.25 * J2, 0.0);
                int d = s ^ (1 << j_0) ^ (1 << j_1);
                int r = 0, l = 0;
                representative(d, &r, &l, N);
                int b = findState(states, r);
                if (b >= 0) {
                    std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N );
                    hamiltonBlock[a][b] += (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * std::exp(numC);
                }
            }
            if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                hamiltonBlock[a][a] += std::complex<double> (0.25 * J2, 0.0);
            } else {
                hamiltonBlock[a][a] -= std::complex<double> (0.25 * J2, 0.0);
                int d = s ^ (1 << j_1) ^ (1 << j_2);
                int r = 0, l = 0;
                representative(d, &r, &l, N);
                int b = findState(states, r);
                if (b >= 0) {
                    std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N );
                    hamiltonBlock[a][b] += (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * std::exp(numC);
                }
            }
        }
    }
}

void momentumBlockSolver(double J1, double J2, int k, const std::vector<int> &states, const std::vector<int> &R_vals, std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXcd> *matrixBlocks) {
    const int statesCount = states.size();
    if (statesCount == 0) {
        return;
    }
    auto **hamiltonBlock = new std::complex<double>*[statesCount];
    for (int i = 0; i < statesCount; i++) {
        hamiltonBlock[i] = new std::complex<double>[statesCount];
        for (int j = 0; j < statesCount; j++) {
            hamiltonBlock[i][j] = 0.0;
        }
    }

    fillHamiltonMomentumBlock(J1, J2, k, states, R_vals, hamiltonBlock);

    Eigen::MatrixXcd H(statesCount, statesCount);
    for (int i = 0; i < statesCount; i++) {
        H.row(i) = Eigen::VectorXcd::Map(&hamiltonBlock[i][0], statesCount);
    }

#if defined(showMatrix) || defined(saveMatrix)
    matrixBlocks->push_back(H);
#endif

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(H);
    const Eigen::VectorXcd& H1EiVal = solver.eigenvalues();
    for (std::complex<double> ev : H1EiVal) {
        HEiValList->push_back(ev);
    }

    for (int i = 0; i < statesCount; i++) {
        delete hamiltonBlock[i];
    }
    delete[] hamiltonBlock;
}

void momentumStateAnsatz(double J1, double J2, std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXcd> *matrixBlocks) {

    int k_lower = -(N+2)/4+1;
    int k_upper = N/4;

    std::vector<std::vector<std::vector<int>>> vec1(N+1, std::vector<std::vector<int>>(N));
    std::vector<std::vector<std::vector<int>>> vec2(N+1, std::vector<std::vector<int>>(N));
    auto *states = &vec1;
    auto *R_vals = &vec2;

    for (int s = 0; s < SIZE; s++) {
        int m = bitSum(s, N);
        for (int k = k_lower; k <= k_upper; k++) {
            int R = checkState(s, k, N);
            if (R >= 0) {
                states->at(m).at(k-k_lower).push_back(s);
                R_vals->at(m).at(k-k_lower).push_back(R);
            }
        }
    }

    for (int m = 0; m <= N; m++) {
        for (int k = k_lower; k <= k_upper; k++) {
            momentumBlockSolver(J1, J2, k, states->at(m).at(k-k_lower), R_vals->at(m).at(k-k_lower), HEiValList, matrixBlocks);
        }
    }

#if defined(showMatrix) || defined(saveMatrix)
    int offset_blocks = 0;
    Eigen::MatrixXcd H_moment_Block = Eigen::MatrixXcd::Zero(SIZE, SIZE);
    for (const Eigen::MatrixXcd& M : *matrixBlocks) {
        H_moment_Block.block(offset_blocks, offset_blocks, M.rows(), M.cols()) = M;
        offset_blocks += M.rows();
    }
#endif
#ifdef showMatrix
    coutMutex.lock();
    std::cout << H_moment_Block << "\n";
    coutMutex.unlock();
#endif
#ifdef saveMatrix
    saveComplexMatrixToFile(H_moment_Block, "HamiltonMomentumStates.txt", "momentum states Ansatz für N = " + std::to_string(N) +
                                                               "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2));
#endif
#if defined(showEigenvalues) || defined(saveEigenvalues)
    HEiValList->shrink_to_fit();
    // sort eigenvalues
    std::sort(HEiValList->begin(), HEiValList->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
        return std::real(c1) < std::real(c2);
    });
#endif
#ifdef showEigenvalues
    coutMutex.lock();
    std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : *HEiValList) {
        std::cout << ev << "\n";
    }
    coutMutex.unlock();
#endif
#ifdef saveEigenvalues
    saveComplexEiVals("EigenvaluesMomentumStates.txt", "momentum states Ansatz für N = " + std::to_string(N) +
                                                      "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *HEiValList);
#endif
}
#endif

/////////////////////////////// parity states (unfinished) ///////////////////////////////
#ifdef parity
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
    saveComplexMatrixToFile(H_parity_Block, "HamiltonParityStates.txt", "parity states Ansatz für N = " + std::to_string(N) +
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
    saveComplexEiVals("EigenvaluesMParityStates.txt", "parity states Ansatz für N = " + std::to_string(N) +
                                                       "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *eiVals);
#endif
}
#endif
/////////////////////////////// MULTITHREADING ///////////////////////////////

#ifdef multiCalc
void threadfunc(double J, int J_pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                double beta, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X) {

    while (true) {

        int p = (int) ( (float) J_pos / (float) J_COUNT * (float) PROGRASSBAR_SEGMENTS);
        coutMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < p; _++) {
            std::cout << "#";
        } for (int _ = p; _ < PROGRASSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) J_pos / (float) J_COUNT * 100.0 ) << "% J1/J2 = " << J << " (" << J_pos << "/" << J_COUNT << ")              ";
        coutMutex.unlock();

        auto *eiVals = new std::vector<std::complex<double>>;
        auto *matrixBlocks = new std::vector<Eigen::MatrixXcd>;
        //auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;

        //naiverAnatz(J, 1.0, eiVals);

        //magnetisierungsAnsatz(J, 1.0, eiVals, matrixBlocks);

        momentumStateAnsatz(J, 1.0, eiVals, matrixBlocks);

        // sort eigenvalues
        eiVals->shrink_to_fit();
        std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });

        // Delta E
        double E0 = std::real(eiVals->at(0));
        double E1 = std::real(eiVals->at(1));

        ///// magnetic susceptibility /////
        double magneticSusceptibility_X;

        // write data
        nextJMutex.lock();
            outDataDeltaE->push_back({J, E1 - E0});
            outDataSpecificHeat_C->push_back({J, getSpecificHeat(beta, *eiVals, N)});
            //outDataMagneticSusceptibility_X->push_back({J, magneticSusceptibility_X});
            J_pos = J_CURRENT;
            J_CURRENT++;
        nextJMutex.unlock();

        if (J_pos > J_COUNT) {
            delete eiVals;
            delete matrixBlocks;
            break;
        } else {
            J = J_START + (J_END-J_START)*J_pos/J_COUNT;
        }
        // clean up
        eiVals->clear();
        matrixBlocks->clear();
    }

}
#endif

/////////////////////////////// MAIN ///////////////////////////////

int main(int argc, char* argv[]) {

    bool silent = false;

    if (argc >= 2) {
        if ( (std::stoi(argv[1]) % 2 == 0) && (std::stoi(argv[1]) >= 6) ) {
             N = std::stoi(argv[1]);
        } else {
            std::cout << "invalid chain size, must be even and at least 6, defaulting to " << N << "\n";
        }
    }
    SIZE = (int) pow(2, N);
    std::cout << "N: " << N << "; size: " << SIZE << std::endl;
    if (argc >= 5) {
        std::cout << "range given: ";
        if (std::stod(argv[2]) > std::stod(argv[3]) || std::stoi(argv[4]) < 1) {
            std::cout << "range invalid, defaulting...\n";
        } else {
            J_START = std::stod(argv[2]); J_END = std::stod(argv[3]); J_COUNT = std::stoi(argv[4]);
        }
        std::cout << "J_START = " << J_START << ", J_END = " << J_END << " and J_COUNT = " << J_COUNT << " from args\n";
    } else {
        std::cout << "no range given: ";
        std::cout << "J_START = " << J_START << ", J_END = " << J_END << " and J_COUNT = " << J_COUNT << " from default\n";
    }

    int cores = cpu_cnt;
    if (argc >= 6) {
        int crs = std::stoi(argv[5]);
        if (crs > 0 && crs <= cpu_cnt) {
            cores = crs;
            std::cout << "using " << cores << "cores\n";
        } else {
        std::cout << "defaulting to using all (" << cores << ") cores\n";
        }
    }

    if (argc >= 7) {
        std::string s1 = "silent";
        std::string s2 = argv[6];
        if (s1 == s2) {
            silent = true;
        }
    }

    // syncing BETA and J ranges
    BETA_START = J_START;
    BETA_END = J_END;
    BETA_COUNT = J_COUNT;

/*
    if (!silent) {
        std::cout << "continue? (y/n):";
        char c;
        std::cin >> c;
        if (c != 'y') {
            wrong_N:
            std::cout << "Enter new N (must be even ans >= 6):";
            int N_usr;
            std::cin >> N_usr;
            if (N_usr >= 6 && N_usr % 2 == 0) {
                N = N_usr;
            } else {
                goto wrong_N;
            }
            wrong_JRANGE:
            std::cout << "Enter new J_START (J1/J2):";
            double JSTART_usr;
            std::cin >> JSTART_usr;
            std::cout << "Enter new J_END (J1/J2):";
            double JEND_usr;
            std::cin >> JEND_usr;
            std::cout << "Enter new J_COUNT (number of datapoints):";
            int JCOUNT_usr;
            std::cin >> JCOUNT_usr;
            if (JSTART_usr <= JEND_usr && JCOUNT_usr >= 1) {
                J_START = JSTART_usr;
                J_END = JEND_usr;
                J_COUNT = JCOUNT_usr;
            } else {
                goto wrong_JRANGE;
            }
        }
    }
    */

/////////////////////////////// MULTITHREADING ///////////////////////////////

#ifdef multiCalc
    const clock_t begin_time = clock();

    std::cout << "\ncalculating:..." << std::endl;

    //auto *outDataDeltaE = new std::vector<std::tuple<double, std::complex<double>>>;
    auto *outDataDeltaE = new std::vector<std::tuple<double, double>>;
    auto *outDataSpecificHeat_C = new std::vector<std::tuple<double, double>>;
    auto *outDataMagneticSusceptibility_X = new std::vector<std::tuple<double, double>>;

    if (J_COUNT < cpu_cnt) {
        cores = J_COUNT;
    }
    std::thread Threads[cores];

    J_CURRENT += cores;

    for (int i = 0; i < cores; i++) {
        Threads[i] = std::thread(threadfunc, J_START + (J_END-J_START)*i/J_COUNT, i + 1, outDataDeltaE, BETA, outDataSpecificHeat_C, outDataMagneticSusceptibility_X);
    }

    for (int i = 0; i < cores; i++) {
        Threads[i].join();
    }

    auto time = float(clock () - begin_time) /  CLOCKS_PER_SEC;
    std::cout << "\n" << "calculations done; this took: " << time << " seconds\n\n";

    // sort datapoints
    std::sort(outDataDeltaE->begin(), outDataDeltaE->end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
        return std::get<0>(a) < std::get<0>(b);
    });
    std::sort(outDataSpecificHeat_C->begin(), outDataSpecificHeat_C->end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
        return std::get<0>(a) < std::get<0>(b);
    });
//    std::sort(outDataMagneticSusceptibility_X->begin(), outDataMagneticSusceptibility_X->end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
//        return std::get<0>(a) < std::get<0>(b);
//    });

    std::string filenameDeltaE = "data_delta_E.txt";
    std::string filenameSpecificHeat_C = "data_specific_heat.txt";
    std::string filenameMagneticSusceptibility_X = "data_magnetic_susceptibility.txt";
    std::string header = "N: " + std::to_string(N) + "\n"
                         + "J1/J2 START: " + std::to_string(J_START) + "\n"
                         + "J1/J2 END: " + std::to_string(J_END) + "\n"
                         + "datapoints: " + std::to_string(J_COUNT) + "\n"
                         + "caculation time with " + std::to_string(cores) + " threads: " + std::to_string(time) + " seconds";

    std::string headerWithBeta = "BETA = " + std::to_string(BETA) + "\n" + header;

    saveOutData(filenameDeltaE, header, "J1/J2", "Delta E in J2", *outDataDeltaE);
    saveOutData(filenameSpecificHeat_C, headerWithBeta, "J1/J2", "specific heat in J2", *outDataSpecificHeat_C);
    //saveOutData(filenameMagneticSusceptibility_X, headerWithBeta, "J1/J2", "magnetic susceptibility in J2", *outDataMagneticSusceptibility_X);

    //return 0;
#endif

/////////////////////////////// naiver Ansatz ///////////////////////////////

#ifdef naiv
    const clock_t begin_time_NAIV = clock();

    std::cout << "\nnaiver Ansatz:..." << std::endl;
    //std::list<std::complex<double>> H_naiv_EiVals;
    auto *H_naiv_EiVals = new std::vector<std::complex<double>>;
    naiverAnatz(J1, J2, H_naiv_EiVals);

    auto time_NAIV = float(clock () - begin_time_NAIV) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_NAIV << " seconds\n";
#endif

/////////////////////////////// fixed magnetization states ///////////////////////////////

#ifdef magnetization
    const clock_t begin_time_MAGNETIZATION = clock();

    std::cout << "\nblockdiagonale m_z:..." << std::endl;
    //std::vector<std::complex<double>> EiVals_m;
    auto *EiVals_m = new std::vector<std::complex<double>>;
    auto *matrixBlocks_m = new std::vector<Eigen::MatrixXd>;
    magnetisierungsAnsatz(J1, J2, EiVals_m, matrixBlocks_m);

    auto time_MAGNETIZATION = float(clock () - begin_time_MAGNETIZATION) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_MAGNETIZATION << " seconds\n";
    delete matrixBlocks_m;
#endif

/////////////////////////////// momentum states ///////////////////////////////

#ifdef momentum
    const clock_t begin_time_MOMENTUM = clock();

    std::cout << "\nmomentum states:..." << std::endl;

    auto *momentEiVals = new std::vector<std::complex<double>>;
    auto *matrixMomentBlocks = new std::vector<Eigen::MatrixXcd>;

    momentumStateAnsatz(J1, J2, momentEiVals, matrixMomentBlocks);

    auto *specificHeat_momentum = new std::vector<std::tuple<double, double>>;

    // specific heat plot
    for (int i = 0; i <= BETA_COUNT; i++) {
        double current_beta = BETA_START + (BETA_END-BETA_START)*i/BETA_COUNT;
        specificHeat_momentum->push_back({current_beta, getSpecificHeat(current_beta, *momentEiVals, N)});
    }

    std::string filenameSpecificHeat_C_momentum = "momentum_specific_heat.txt";
    std::string header_momentum = "N: " + std::to_string(N) + "\n"
                                + "BETA START: " + std::to_string(BETA_START) + "\n"
                                + "BETA END: " + std::to_string(BETA_END) + "\n"
                                + "datapoints: " + std::to_string(BETA_COUNT) + "\n"
                                + "caculation time with " + std::to_string(cores) + " threads: " + std::to_string(time) + " seconds";

    std::string headerWithJ_momentum = "J1/J2 = " + std::to_string(J1/J2) +"\n" + header;
    saveOutData(filenameSpecificHeat_C_momentum, headerWithJ_momentum, "J1/J2", "specific heat in J2", *specificHeat_momentum);

    auto time_MOMENTUM = float(clock () - begin_time_MOMENTUM) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_MOMENTUM << " seconds\n";
    delete momentEiVals;
    delete matrixMomentBlocks;
#endif

/////////////////////////////// parity states (unfinished) ///////////////////////////////

#ifdef parity
    const clock_t begin_time_PARITY = clock();

    std::cout << "\nparity states:..." << std::endl;

    auto *parityEiVals = new std::vector<std::complex<double>>;
    auto *matrixParityBlocks = new std::vector<Eigen::MatrixXcd>;

    parityStateAnsatz(J1, J2, parityEiVals, matrixParityBlocks);

    auto time_PARITY = float(clock () - begin_time_PARITY) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_PARITY << " seconds\n";
    delete parityEiVals;
    delete matrixParityBlocks;
#endif

    return 0;
}
