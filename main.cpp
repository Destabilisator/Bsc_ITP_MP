#include "main.h"

//// methods ////
#define naiv
#define magnetization
#define momentum

#define multiCalc

///// output /////
//#define showMatrix
//#define saveMatrix
//#define showEigenvalues
//#define saveEigenvalues

///// global variables /////
int N = 10; // has to be even to preserve the periodic boundary conditions of the delta chain
int SIZE;
double J_START, J_END;
int J_COUNT;
int J_CURRENT = 0;

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

void naiverAnatz(double J1, double J2, std::list<std::complex<double>> &HEiValList) {
    std::cout << "\nnaiver Ansatz:..." << std::endl;
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
    std::cout << "solving...\n";
    Eigen::EigenSolver<Eigen::MatrixXd> solver1(H1);
    const Eigen::VectorXcd &H1EiVal = solver1.eigenvalues();
    for (std::complex<double> ev : H1EiVal) {
        HEiValList.push_back(ev);
    }
    // sort List
    HEiValList.sort([](const std::complex<double> &c1, const std::complex<double> &c2) {
       return std::real(c1) < std::real(c2);
    });
#ifdef showEigenvalues
    std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : HEiValList) {
        std::cout << ev << "\n";
    }
#endif
#ifdef saveEigenvalues
    saveComplexEiVals("EigenvaluesNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), HEiValList);
#endif
    for (int i = 0; i < SIZE; i++) {
        delete hamilton1[i];
    } delete[] hamilton1;
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
                //continue;
                int d = s ^ (1 << j_0) ^ (1 << j_2);
                int pos_d = findState(states, d);
                hamiltonBlock[i][pos_d] = 0.5 * J1;
            }
            if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                hamiltonBlock[i][i] += 0.25 * J2;
            } else {
                hamiltonBlock[i][i] -= 0.25 * J2;
                //continue;
                int d = s ^ (1 << j_0) ^ (1 << j_1);
                int pos_d = findState(states, d);
                hamiltonBlock[i][pos_d] = 0.5 * J2;
            }
            if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                hamiltonBlock[i][i] += 0.25 * J2;
            } else {
                hamiltonBlock[i][i] -= 0.25 * J2;
                //continue;
                int d = s ^ (1 << j_1) ^ (1 << j_2);
                int pos_d = findState(states, d);
                hamiltonBlock[i][pos_d] = 0.5 * J2;
            }
        }
    }
}

void magBlock_getEiVal(double J1, double J2, int m, std::vector<std::complex<double>> &HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks) {
//    coutMutex.lock();
//    std::cout << "number of up spins: " << m << " out of " << N << "\n";
//    coutMutex.unlock();
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
    //EiValWriterMutex.lock();
    for (std::complex<double> ev : H1EiVal) {
        HEiValList.push_back(ev);
    }
    //EiValWriterMutex.unlock();

    states->clear();
    delete states;
    for (int i = 0; i < statesCount; i++) {
        delete hamiltonBlock[i];
    }
    delete[] hamiltonBlock;
}

void magnetisierungsAnsatz(double J1, double J2, std::vector<std::complex<double>> &HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks) {
    for (int m = 0; m <= N; m++) {
        magBlock_getEiVal(J1, J2, m, HEiValList, matrixBlocks);
    }
    // sort eigenvalues
    std::sort(HEiValList.begin(), HEiValList.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
        return std::real(c1) < std::real(c2);
    });
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
#ifdef showEigenvalues
    coutMutex.lock();
    std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : HEiValList) {
        std::cout << ev << "\n";
    }
    coutMutex.unlock();
#endif
#ifdef saveEigenvalues
    saveComplexEiVals("EigenvaluesMagnetization.txt", "Magnetisierungs Ansatz für N = " + std::to_string(N) +
                      "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), HEiValList);
#endif
}
#endif

/////////////////////////////// momentum states (unfinished) ///////////////////////////////

#ifdef momentum

void fillHamiltonMomentumBlock(double J1, double J2, int k,const std::vector<int> &states, const std::vector<int> &R_vals, std::complex<double> **hamiltonBlock) {
    //std::cout << "states.size(): " << states.size() << "\n";
    for (int a = 0; a < states.size(); a++) {
        int s = states.at(a);
        std::cout << "state: ";
        printBits(s, N);
        for (int n = 0; n < N/2; n++) {
            //std::cout << "Periodizitaeten: Ra " << R_vals.at(a);
            // declaring indices
            int j_0 = 2 * n;
            int j_1 = (j_0+1) % N;
            int j_2 = (j_0+2) % N;
            // applying H to state s
            std::cout << "j_0 = " << j_0 << ", j_2 = " << j_2 << "\n";
            if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
//                std::cout << "diagonal ";
//                printBits(s, N);
                std::cout << "changes H(" << a << ", " << a << ") from " << hamiltonBlock[a][a];
                hamiltonBlock[a][a] += std::complex<double> (0.25 * J1, 0.0);
                std::cout << " to " << hamiltonBlock[a][a] << "\n";
            } else {
                std::cout << "changes H(" << a << ", " << a << ") from " << hamiltonBlock[a][a];
                hamiltonBlock[a][a] -= std::complex<double> (0.25 * J1, 0.0);
                std::cout << " to " << hamiltonBlock[a][a] << "\n";
                int d = s ^ (1 << j_0) ^ (1 << j_2);
                int r = 0, l = 0;
                representative(d, &r, &l, N);
                int b = findState(states, r);
                if (b >= 0) {
                    std::cout << "off diagonal: ";
                    printBits(r, N);
                    //std::cout << ", Rb " << R_vals.at(b);
                    //continue;
                    //std::complex<double> numC(0.0, PI * (double) k * (double) l / (double) N );
                    std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N );
                    std::cout << "changes H(" << a << ", " << b << ") from " << hamiltonBlock[a][b];
                    hamiltonBlock[a][b] += (std::complex<double>) 0.5 * J1 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * std::exp(numC);
                    std::cout << " to " << hamiltonBlock[a][b] << "\n";
                    //std::cout << (std::complex<double>) 0.5 * J1 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * exp(numC) << "\n";
                }
            }
            std::cout << "j_0 = " << j_0 << ", j_1 = " << j_1 << "\n";
            if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                std::cout << "changes H(" << a << ", " << a << ") from " << hamiltonBlock[a][a];
                hamiltonBlock[a][a] += std::complex<double> (0.25 * J2, 0.0);
                std::cout << " to " << hamiltonBlock[a][a] << "\n";
            } else {
                std::cout << "changes H(" << a << ", " << a << ") from " << hamiltonBlock[a][a];
                hamiltonBlock[a][a] -= std::complex<double> (0.25 * J2, 0.0);
                std::cout << " to " << hamiltonBlock[a][a] << "\n";
                int d = s ^ (1 << j_0) ^ (1 << j_1);
                int r = 0, l = 0;
                representative(d, &r, &l, N);
                int b = findState(states, r);
                if (b >= 0) {
                    std::cout << "off diagonal: ";
                    printBits(r, N);
                    //std::cout << ", Rb " << R_vals.at(b);
                    //continue;
                    //std::complex<double> numC(0.0, PI * (double) k * (double) l / (double) N );
                    std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N );
                    std::cout << "changes H(" << a << ", " << b << ") from " << hamiltonBlock[a][b];
                    hamiltonBlock[a][b] += (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * std::exp(numC);
                    std::cout << " to " << hamiltonBlock[a][b] << "\n";
                    //std::cout << (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * exp(numC) << "\n";
                }
            }
            std::cout << "j_1 = " << j_1 << ", j_2 = " << j_2 << "\n";
            if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                std::cout << "changes H(" << a << ", " << a << ") from " << hamiltonBlock[a][a];
                std::cout << " to " << hamiltonBlock[a][a] << "\n";
                hamiltonBlock[a][a] += std::complex<double> (0.25 * J2, 0.0);
            } else {
                std::cout << "changes H(" << a << ", " << a << ") from " << hamiltonBlock[a][a];
                hamiltonBlock[a][a] -= std::complex<double> (0.25 * J2, 0.0);
                std::cout << " to " << hamiltonBlock[a][a] << "\n";
                int d = s ^ (1 << j_1) ^ (1 << j_2);
                int r = 0, l = 0;
                representative(d, &r, &l, N);
                int b = findState(states, r);
                if (b >= 0) {
                    std::cout << "off diagonal: ";
                    printBits(r, N);
                    //std::cout << ", Rb " << R_vals.at(b);
                    //continue;
                    //std::complex<double> numC(0.0, PI * (double) k * (double) l / (double) N );
                    std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N );
                    std::cout << "changes H(" << a << ", " << b << ") from " << hamiltonBlock[a][b];
                    hamiltonBlock[a][b] += (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * std::exp(numC);
                    std::cout << " to " << hamiltonBlock[a][b] << "\n";
                    //std::cout << (std::complex<double>) 0.5 * J2 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * exp(numC) << "\n";
                }
            }
            //std::cout << "\n";
        }
    }
}

void momentumBlockSolver(double J1, double J2, int k, const std::vector<int> &states, const std::vector<int> &R_vals, std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXcd> *matrixBlocks) {
    const int statesCount = states.size();
    //std::cout << "statesCount: " << statesCount << " " << R_vals.size() << "\n";
    if (statesCount == 0) {
        //std::cout << "empty block\n";
        return;
    }
    auto **hamiltonBlock = new std::complex<double>*[statesCount];
    for (int i = 0; i < statesCount; i++) {
        hamiltonBlock[i] = new std::complex<double>[statesCount];
        for (int j = 0; j < statesCount; j++) {
            hamiltonBlock[i][j] = 0.0;
        }
    }
    //std::cout << "filling block\n";
    fillHamiltonMomentumBlock(J1, J2, k, states, R_vals, hamiltonBlock);

    Eigen::MatrixXcd H(statesCount, statesCount);
    for (int i = 0; i < statesCount; i++) {
        H.row(i) = Eigen::VectorXcd::Map(&hamiltonBlock[i][0], statesCount);
    }

#if defined(showMatrix) || defined(saveMatrix)
    matrixBlocks->push_back(H);
#endif

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(H);
    //Eigen::EigenSolver<Eigen::MatrixXcd> solver(H);
    const Eigen::VectorXcd& H1EiVal = solver.eigenvalues();
    //EiValWriterMutex.lock();
    for (std::complex<double> ev : H1EiVal) {
        HEiValList->push_back(ev);
    }
    //EiValWriterMutex.unlock();

    for (int i = 0; i < statesCount; i++) {
        delete hamiltonBlock[i];
    }
    delete[] hamiltonBlock;
}

void momentumStateAnsatz(double J1, double J2, std::vector<std::complex<double>> *HEiValList, std::vector<Eigen::MatrixXcd> *matrixBlocks) {
    std::vector<std::vector<std::vector<int>>> vec1(N+1, std::vector<std::vector<int>>(N));
    std::vector<std::vector<std::vector<int>>> vec2(N+1, std::vector<std::vector<int>>(N));
    auto *states = &vec1;
    auto *R_vals = &vec2;
    for (int s = 0; s < SIZE; s++) {
        int m = bitSum(s, N);
        for (int k = -N/4 + 1; k <= N/4; k++) { // for (int k = -N/2 + 1; k <= N/2; k++) {
            int R = checkState(s, k, N);
            if (R >= 0) {
//                std::cout << "m = " << m << " k = " << k << "\n";
//                std::cout << "R = " << R  << ": ";
//                printBits(s, N);
                states->at(m).at(k+N/4-1).push_back(s); // states->at(m).at(k+N/2-1).push_back(s);
                R_vals->at(m).at(k+N/4-1).push_back(R); // R_vals->at(m).at(k+N/2-1).push_back(R);
            }
        }
    }
//        for (int k = -N/2 + 1; k <= N/2; k++) { // for (int k = -N/2 + 1; k <= N/2; k++) {
//            int R = checkState(s, k, N);
//            if (R >= 0) {
////                std::cout << "m = " << m << " k = " << k << "\n";
////                std::cout << "R = " << R  << ": ";
////                printBits(s, N);
//                states->at(m).at(k+N/2-1).push_back(s); // states->at(m).at(k+N/2-1).push_back(s);
//                R_vals->at(m).at(k+N/2-1).push_back(R); // R_vals->at(m).at(k+N/2-1).push_back(R);
//            }
//        }
//    }


    std::cout << "calculating eigenvalues\n";
    for (int m = 0; m <= N; m++) {
        for (int k = -N/4 + 1; k <= N/4; k++) {
        //for (int k = -N/2 + 1; k <= N/2; k++) {
            std::cout << "m = " << m << " k = " << k << "\n";
            momentumBlockSolver(J1, J2, k, states->at(m).at(k+N/4-1), R_vals->at(m).at(k+N/4-1), HEiValList, matrixBlocks);
            //momentumBlockSolver(J1, J2, k, states->at(m).at(k+N/2-1), R_vals->at(m).at(k+N/2-1), HEiValList, matrixBlocks);
        }
    }

    // sort eigenvalues
    std::sort(HEiValList->begin(), HEiValList->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
        return std::real(c1) < std::real(c2);
    });

#if defined(showMatrix) || defined(saveMatrix)
    int offset_mag_blocks = 0;
    Eigen::MatrixXcd H_moment_Block = Eigen::MatrixXcd::Zero(SIZE, SIZE);
    for (const Eigen::MatrixXcd& M : *matrixBlocks) {
        H_moment_Block.block(offset_mag_blocks, offset_mag_blocks, M.rows(), M.cols()) = M;
        offset_mag_blocks += M.rows();
    }
#endif
#ifdef showMatrix
    coutMutex.lock();
    std::cout << H_moment_Block << "\n";
    coutMutex.unlock();
#endif
#ifdef saveMatrix
    saveComplexMatrixToFile(H_moment_Block, "HamiltonMomentumStates.txt", "momentum state Ansatz für N = " + std::to_string(N) +
                                                               "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2));
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
    saveComplexEiVals("EigenvaluesMomentumStates.txt", "momentum state Ansatz für N = " + std::to_string(N) +
                                                      "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *HEiValList);
#endif
}

#endif

/////////////////////////////// MULTITHREADING ///////////////////////////////

#ifdef multiCalc
void threadfunc(double J, std::vector<std::tuple<double, std::complex<double>>> *outData, int J_pos) {

    coutMutex.lock();
    std::cout << "J1/J2 = " << J << " (" << J_pos << "/" << J_COUNT << ")\n";
    coutMutex.unlock();

    //auto *eiVals = new std::vector<std::complex<double>>;
    //auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;
    std::vector<std::complex<double>> eiVals;
    auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;

    magnetisierungsAnsatz(J, 1, eiVals, matrixBlocks);
    // sort eigenvalues
    std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
        return std::real(c1) < std::real(c2);
    });

    std::complex<double> E0 = eiVals.at(0);
    std::complex<double> E1;

    for (int i = 1; i < eiVals.size(); i++) {
        if (abs(E0 - eiVals.at(i)) > 0.001) {
            E1 = eiVals.at(i);
            break;
        }
    }

    nextJMutex.lock();
    outData->push_back({J, E1 - E0});
    J_pos = J_CURRENT;
    J_CURRENT++;
    nextJMutex.unlock();

    if (J_pos > J_COUNT) {
        return;
    } else {
        threadfunc( J_START + (J_END-J_START)*J_pos/J_COUNT, outData, J_pos);
    }
}
#endif

/////////////////////////////// MAIN ///////////////////////////////

int main(int argc, char* argv[]) {

    if (argc >= 2) {
        std::cout << "N from args\n";
        if ( (std::stoi(argv[1]) % 2 == 0) && (std::stoi(argv[1]) >= 4) ) {
             N = std::stoi(argv[1]);
        } else {
            std::cout << "invalid chain SIZE, must be even and at least 4, defaulting to " << N << "\n";
        }
    }
    if (argc == 5) {
        std::cout << "range given: ";
        J_START = std::stod(argv[2]); J_END = std::stod(argv[3]); J_COUNT = std::stoi(argv[4]);
        if (J_START > J_END | J_COUNT < 1) {
            std::cout << "range invalid, defaulting...\n";
            goto default_J;
        }
        std::cout << "J_START = " << J_START << ", J_END = " << J_END << " and J_COUNT = " << J_COUNT << " from args\n";
    } else {
        default_J:
        J_START = 0.8; J_END = 1.4; J_COUNT = 2000;
        std::cout << "J_START = " << J_START << ", J_END = " << J_END << " and J_COUNT = " << J_COUNT << " from default\n";
    }

    SIZE = (int) pow(2, N);
    std::cout << "N: " << N << "; SIZE: " << SIZE << std::endl;

/////////////////////////////// MULTITHREADING ///////////////////////////////

#ifdef multiCalc
    const clock_t begin_time = clock();

    std::cout << "\ncalculating:..." << std::endl;

    auto *outData = new std::vector<std::tuple<double, std::complex<double>>>;

    int cores = cpu_cnt;
    if (J_COUNT < cpu_cnt) {
        cores = J_COUNT;
    }
    std::thread Threads[cores];

    J_CURRENT += cores;

    for (int i = 0; i < cores; i++) {
//        coutMutex.lock();
//        std::cout << "starting thread " << i << "\n";
//        coutMutex.unlock();
        Threads[i] = std::thread(threadfunc, J_START + (J_END-J_START)*i/J_COUNT, outData, i+1);
    }

    for (int i = 0; i < cores; i++) {
        Threads[i].join();
    }

    auto time = float(clock () - begin_time) /  CLOCKS_PER_SEC;
    std::cout << "\n" << "calculations done; this took: " << time << " seconds\n";

    std::string filename = "data.txt";//_" + std::to_string(N) + "_" + std::to_string(J_START) + "_" + std::to_string(J_END) + "_" + std::to_string(J_COUNT) + ".txt";
    std::cout << "saving to file " << filename << "...\n";

    // sort datapoints
    std::sort(outData->begin(), outData->end(), [](const std::tuple<double, std::complex<double>> &a, const std::tuple<double, std::complex<double>> &b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    std::ofstream file;
    try {
        file.open("./results/" + filename);
        file << "N: " << N << "\n";
        file << "J1/J2 START_ " << J_START << "\n";
        file << "J1/J2 END: " << J_END << "\n";
        file << "datapoints: " << J_COUNT << "\n";
        file << "caculation time with " << cores << " threads: " << time << "\n\n";
        file << "J1/J2\tDelta E in J2\n";
        for (std::tuple<double, std::complex<double>> data : *outData) {
            file << std::get<0>(data) << "\t" << std::abs(std::get<1>(data)) << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
        return 1;
    }
    file.close();

    return 0;
#endif

    const double J1 = 1.0, J2 = 0.0;

/////////////////////////////// naiver Ansatz ///////////////////////////////

#ifdef naiv
    const clock_t begin_time_NAIV = clock();

    std::list<std::complex<double>> H_naiv_EiVals;
    naiverAnatz(J1, J2, H_naiv_EiVals);

    auto time_NAIV = float(clock () - begin_time_NAIV) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_NAIV << " seconds\n";
#endif

/////////////////////////////// fixed magnetization states ///////////////////////////////

#ifdef magnetization
    const clock_t begin_time_MAGNETIZATION = clock();

    std::cout << "\nblockdiagonale m_z:..." << std::endl;
    std::vector<std::complex<double>> H_mag_EiVals;
    auto *matrixBlocks_m = new std::vector<Eigen::MatrixXd>;
    magnetisierungsAnsatz(J1, J2, H_mag_EiVals, matrixBlocks_m);

    auto time_MAGNETIZATION = float(clock () - begin_time_MAGNETIZATION) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_MAGNETIZATION << " seconds\n";
    delete matrixBlocks_m;
#endif

/////////////////////////////// momentum states (unfinished) ///////////////////////////////

#ifdef momentum
    const clock_t begin_time_MOMENTUM = clock();

    std::cout << "\nmomentum states:..." << std::endl;

    auto *H_moment_EiVals = new std::vector<std::complex<double>>;
    auto *matrixBlocks = new std::vector<Eigen::MatrixXcd>;

    momentumStateAnsatz(J1, J2, H_moment_EiVals, matrixBlocks);

    auto time_MOMENTUM = float(clock () - begin_time_MOMENTUM) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_MOMENTUM << " seconds\n";
    delete matrixBlocks;
#endif

    return 0;
}
