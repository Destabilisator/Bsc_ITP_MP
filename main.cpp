#include "main.h"

//// methods ////
//#define naiv
#define magnetization
//#define momentum

//// output ////
//#define showMatrix
//#define saveMatrix
#define showEigenvalues
//#define saveEigenvalues

//// global variables ////
int N = 2*2; // has to be even to preserve the periodic boundary conditions of the delta chain
int size;
double J_start, J_end;
int J_count, J_current;

/////////////////////////////// naiver Ansatz ///////////////////////////////

#ifdef naiv
void fillHamiltonNaiv(double** hamilton, double J1, double J2) {
    for (int s = 0; s <= size -1; s++) {
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
    static auto **hamilton1 = new double*[size];
    for (int i = 0; i < size; i++) {
        hamilton1[i] = new double[size];
        for (int j = 0; j < size; j++) {
            hamilton1[i][j] = 0.0;
        }
    }
    fillHamiltonNaiv(hamilton1, J1, J2);
    Eigen::MatrixXd H1(size, size);
    for (int i = 0; i < size; i++) {
        H1.row(i) = Eigen::VectorXd::Map(&hamilton1[i][0], size);
    }
#ifdef showMatrix
    std::cout << H1 << std::endl;
#endif
#ifdef saveMatrix
    saveHamilton(hamilton1, "HamiltonNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), size);
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
    for (int i = 0; i < size; i++) {
        delete hamilton1[i];
    } delete[] hamilton1;
}
#endif

/////////////////////////////// fixed magnetization states ///////////////////////////////

#ifdef magnetization
void fillHamiltonBlock(double J1, double J2, const std::vector<int>& states, double **hamiltonBlock) {
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

void magBlock_getEiVal(double J1, double J2, int m, std::list<std::complex<double>> &HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks) {
//    coutMutex.lock();
//    std::cout << "number of up spins: " << m << " out of " << N << "\n";
//    coutMutex.unlock();
    auto *states = new std::vector<int>;
    fillStates(states, m, N, size);
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

void magnetisierungsAnsatz(double J1, double J2, std::list<std::complex<double>> &HEiValList, std::vector<Eigen::MatrixXd> *matrixBlocks) {
    for (int m = 0; m <= N; m++) {
        magBlock_getEiVal(J1, J2, m, HEiValList, matrixBlocks);
    }
    // sort List
    HEiValList.sort([](const std::complex<double> &c1, const std::complex<double> &c2) {
        return std::real(c1) < std::real(c2);
    });
#if defined(showMatrix) || defined(saveMatrix)
    int offset_mag_blocks = 0;
    Eigen::MatrixXd H_mag_Block = Eigen::MatrixXd::Zero(size, size);
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
                      "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), H_mag_EiVals);
#endif
}
#endif

/////////////////////////////// momentum states (unfinished) ///////////////////////////////

#ifdef momentum


#endif

/////////////////////////////// MAIN ///////////////////////////////

int main(int argc, char* argv[]) {

    int J1 = 1, J2 = 1;

    if (argc >= 2) {
        std::cout << "N from args\n";
        if ( (std::stoi(argv[1]) % 2 == 0) && (std::stoi(argv[1]) >= 4) ) {
             N = std::stoi(argv[1]);
        } else {
            std::cout << "invalid chain size, must be even and at least 4, defaulting to " << N << "\n";
        }
    }
    if (argc == 5) {
        std::cout << "range given: ";
        J_start = std::stod(argv[2]); J_end = std::stod(argv[3]); J_count = std::stoi(argv[4]);
        if (J_start > J_end | J_count < 1) {
            std::cout << "range invalid, defaulting...\n";
            goto default_J;
        }
        std::cout << "J_start = " << J_start << ", J_end = " << J_end << " and J_count = " << J_count << " from args\n";
    } else {
        default_J:
        J_start = 1; J_end = 1; J_count = 1;
        std::cout << "J_start = " << J_start << ", J_end = " << J_end << " and J_count = " << J_count << " from default\n";
    }

    size = (int) pow(2, N);
    std::cout << "N: " << N << "; size: " << size << std::endl;

/////////////////////////////// naiver Ansatz ///////////////////////////////

#ifdef naiv
    const clock_t begin_time_NAIV = clock();

    std::list<std::complex<double>> H_naiv_EiVals;
    naiverAnatz(J1, J2, H_naiv_EiVals);

    auto time_NAIV = float( clock () - begin_time_NAIV ) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_NAIV << " seconds\n";
#endif

/////////////////////////////// fixed magnetization states ///////////////////////////////

#ifdef magnetization
    const clock_t begin_time_MAGNETIZATION = clock();

    std::cout << "\nblockdiagonale m_z:..." << std::endl;
    std::list<std::complex<double>> H_mag_EiVals;
    auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;
    magnetisierungsAnsatz(J1, J2, H_mag_EiVals, matrixBlocks);

    auto time_MAGNETIZATION = float( clock () - begin_time_MAGNETIZATION ) /  CLOCKS_PER_SEC;
    std::cout << "calculations done; this took: " << time_MAGNETIZATION << " seconds\n";
    delete matrixBlocks;
#endif

/////////////////////////////// momentum states (unfinished) ///////////////////////////////

#ifdef momentum
    // Using momentum states.
    std::cout << "\nmomentum states:..." << std::endl;

    //std::cout << "allocating matrix\n";
    static auto **hamilton3 = new std::complex<double>*[size];
    for (int i = 0; i < size; i++) {
        hamilton3[i] = new std::complex<double>[size];
        for (int j = 0; j < size; j++) {
            hamilton3[i][j] = std::complex<double> (0.0, 0.0);
        }
    }


    auto *statesM = new std::vector<int>;
    auto *statesListM = new std::vector<int>;
    auto *statesPerioM = new std::vector<int>;
    offset = 0;
    for (int k = -(N-1)/2; k <= N/2 ; k++) {
        //std::cout << "k: " << k << "\n";
        for (int m = 0; m <= N; m++) {
            fillStates(statesM, m);
            for (int a : *statesM) {
                //std::cout << "found state: ";
                //printBits(a);
                int perio = checkState(a, k);
                if (perio >= 0) {
                    //std::cout << perio << "\n";
                    statesListM->push_back(a);
                    //std::cout << statesList->at(0) << ": ";
                    statesPerioM->push_back(perio);
                    //printBits(a);
                    //std::cout << statesPerio->at(0) << "\n";
                }
            }
            int M = statesListM->size();
            auto **hamiltonBlock = new std::complex<double>*[M];
            for (int i = 0; i < M; i++) {
                hamiltonBlock[i] = new std::complex<double>[M];
                for (int j = 0; j < M; j++) {
                    hamiltonBlock[i][j] = std::complex<double> (0.0, 0.0);
                }
            }
            //std::cout << "filling hamilton block\n";
            fillHamiltonMomentumBlock(hamiltonBlock, statesListM, statesPerioM, k);
            //std::cout << "writing hamilton block to full\n";
            writeHamiltonMomentumBlockToFull(hamiltonBlock, hamilton3, M, offset);
            offset += M;

            //std::cout << "clean up\n";
            statesM->clear();
            statesListM->clear();
            statesPerioM->clear();
            for (int i = 0; i < M; i++) {
                delete hamiltonBlock[i];
            } delete[] hamiltonBlock;
        }
    }

    Eigen::MatrixXcd H3(size, size);
    for (int i = 0; i < size; i++) {
        H3.row(i) = Eigen::RowVectorXcd::Map(&hamilton3[i][0], size);
    }
#ifdef showMatrix
    std::cout << H3 << std::endl;
#endif
#ifdef calculateEigenvalues
    std::cout << "solving...\n";
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver3(H3);
    Eigen::VectorXcd H3EiVal = solver3.eigenvalues();
    std::cout << H3EiVal << std::endl;
    saveComplexHamilton(hamilton3, "Hamilton3.txt", "momentum states", H3EiVal);
#endif
    for (int i = 0; i < size; i++) {
        delete hamilton3[i];
    } delete[] hamilton3;
#endif


    return 0;
}
