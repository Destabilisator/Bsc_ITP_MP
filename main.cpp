//#include "main.h"

#include <iostream>
#include <fstream>
#include <bitset>
#include <cmath>
#include <string>
#include <vector>
#include <Eigen/Eigenvalues>
#include <iterator>
#include <list>
#include "headers/helpers.h"
//#include <algorithm>

#define PI  3.14159265358979323846

// output
#define calculateEigenvalues
#define saveOnlyEiVal
#define showMatrix

// methods
#define naiv
#define magnetization
//#define momentum

int N = 2*2; // has to be odd to conserve delta shape
int size;

/////////////////////////////// naiver Ansatz ///////////////////////////////
// naiver Ansatz: direkt füllen, ohne zu ordnen
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
    saveHamilton(hamilton1, "HamiltonNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), size);
#endif
#ifdef calculateEigenvalues
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
    saveComplexEiVals("EigenvaluesNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), HEiValList);
#endif
    for (int i = 0; i < size; i++) {
        delete hamilton1[i];
    } delete[] hamilton1;
}

void magBlock_getEiVal(double J1, double J2, int m, std::list<std::complex<double>> &HEiValList) {

}

/////////////////////////////// MAIN ///////////////////////////////
int main(int argc, char* argv[]) {

    if (argc == 2) {
        std::cout << "N from args\n";
        if (std::stoi(argv[1]) % 2 == 0 | std::stoi(argv[1]) < 4) {
             N = std::stoi(argv[1]);
        } std::cout << "invalid chain size, must be even... defaulting to " << N << "\n";
    }
    std::cout << "no or invalid args given, default N\n";
    size = (int) pow(2, N);
    std::cout << "N: " << N << "; size: " << size << std::endl;

/////////////////////////////// naiver Ansatz ///////////////////////////////

#ifdef naiv
    std::list<std::complex<double>> H_naiv_EiVals;
    naiverAnatz(1, 1, H_naiv_EiVals);
    std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : H_naiv_EiVals) {
        std::cout << ev << "\n";
    }
#endif

    //int offset;

/////////////////////////////// fixed magnetization states (unfinished) ///////////////////////////////
#ifdef magnetization
    std::cout << "\nblockdiagonale m_z" << std::endl;
    static auto **hamilton2 = new float*[size];
    for (int i = 0; i < size; i++) {
        hamilton2[i] = new float[size];
        for (int j = 0; j < size; j++) {
            hamilton2[i][j] = 0.0;
        }
    }

    //std::vector<int> statesBlocks[N+1];
    //offset = 0;
    std::list<std::complex<double>> H_mag_EiVals;

    for (int m = 0; m <= N; m++) {
#ifdef calculateEigenvalues
        magBlock_getEiVal(1, 1, m, H_mag_EiVals);
#endif
    }

    Eigen::MatrixXf H2(size, size);
    for (int i = 0; i < size; i++) {
        H2.row(i) = Eigen::VectorXf::Map(&hamilton2[i][0], size);
    }
#ifdef showMatrix
    std::cout << H2 << std::endl;
#endif
#ifdef calculateEigenvalues
//    std::cout << "solving...\n";
//    Eigen::EigenSolver<Eigen::MatrixXf> solver2(H2);
//    Eigen::VectorXcd H2EiVal = solver2.eigenvalues();
//    std::cout << H2EiVal << std::endl;
//    saveHamilton(hamilton2, "Hamilton2.txt", "Blöcke konstanter Magnetisierung für N = " + std::to_string(N), H2EiVal);
#endif
    for (int i = 0; i < size; i++) {
        delete hamilton2[i];
    } delete[] hamilton2;
#endif


/////////////////////////////// momentum states (not started) ///////////////////////////////
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
