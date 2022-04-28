#include <iostream>
#include <fstream>
#include <bitset>
#include <math.h>
#include <string>
#include <vector>
#include <Eigen/Eigenvalues>
#include <iterator>
#include <list>

#define PI  3.14159265358979323846

// output
#define calculateEigenvalues
#define saveOnlyEiVal
#define showMatrix

// methods
#define naiv
//#define magnetization
//#define momentum

int N = 3;
int size;

void printBits(int a) {
    for (int n = N-1; n >= 0; n--) {
        std::cout << ( (a >> n) & 1 );
    } std::cout << "\n";
}

// cyclicly translate bits in s by n to the right
int translateRight(int s, int n) {
    for (int _ = 0; _ < n; _++) {
        int bit = s & 1;
        s = (s>>1) | (bit << (N-1));
    } return s;
}

// cyclicly translate bits in s by n to the left
int translateLeft(int s, int n) {
    for (int _ = 0; _ < n; _++) {
        int bit = (s >> (N-1)) & 1;
        s = (s<<1) | bit;
        s &= ~(1 << N);
    } return s;
}

// spiegelt die bits
int reflectBits(int s) {
    int t = 0;
    for (int i = 0; i < N; i++) {
        t |= ((s >> (N - 1 - i)) & 1) << i;
    }
    return t;
}

int invertBits(int s) {
    return s ^ (int) pow(2,N)-1;
}

// sum up all bits in s
int bitSum(int s) {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        sum += ((s >> i) & 1 );
    } return sum;
}

void saveHamilton(double** hamilton, const std::string &filename, const std::string &header) {
    std::cout << "saving to file '" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open(filename);
        file << header << "\n";
        for (int i = 0; i <= size -1; i++) {
            for (int j = 0; j <= size-1; j++) {
                if (hamilton[i][j] < 0.001 && hamilton[i][j] > -0.001) {
                    file << " \t";
                } else {
                    file << hamilton[i][j] << "\t" ;
                }
            }
            file << std::endl;
        }

    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveEiVals(const std::string &filename, const std::string &header, const std::list<double> &eiVals) {
    std::cout << "saving to file '" << filename << "'..." << std::endl;
    std::ofstream file;
    try {
        file.open(filename);
        file << header << "\n";
        file << "\nEigenvalues:\n";
        for (double ev : eiVals) {
            file << ev << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    std::cout << "done\n";
    file.close();
}

// naiver Ansatz: direkt füllen, ohne zu ordnen
void fillHamiltonNaiv(double** hamilton, double J1, double J2) {
    for (int s = 0; s <= size -1; s++) {
        for (int i = 0; i <= N-1; i++) {
            int j = (i+2) % N;
            if (((s >> i) & 1) == ((s >> j) & 1)) {
                hamilton[s][s] += 0.25 * J1;
            } else {
                hamilton[s][s] -= 0.25 * J1;
                int b = s ^ (1 << i) ^ (1 << j);
                //std::cout << "off " <<  a << " " << b << "\n";
                hamilton[s][b] = 0.5;
            }
        }
    }
}

// findet alle Zustände mit der Magnetisierung m_z = m
void fillStates(std::vector<int> *states, int m) {
    for (int s = 0; s <= size - 1; s++) {
        if (bitSum(s) == m) {
            //std::cout << bitSum(s) << ": ";
            states->push_back(s);
            //printBits(s);
        }
    } states->shrink_to_fit();
}

// findet einen Zustand in der Zustandsliste
int findState(std::vector<int> *states, int s, int M) {
    int pos, pos_min = 0, pos_max = M-1;
    while (true) {
        pos = pos_min + (pos_max - pos_min ) / 2;
        if (s < states->at(pos)) {
            pos_max = pos - 1;
        } else if (s > states->at(pos)) {
            pos_min = pos + 1;
        } else {
            return pos;
        } if (pos_min > pos_max) {
            return -1;
        }
    }
}

// füllt einen Hamiltonblock aus Zuständen mit gegebener Magnetisierung aus
void fillHamiltonBlock(std::vector<int> *states, float** hamilton, int M) {
    for (int k = 0; k < M; k++) {
        //std::cout << k << ": ";
        int a = states->at(k);
        //printBits(a);
        for (int i = 0; i <= N-1; i++) {
            int j = (i+1) % N;
            if (((a >> i) & 1) == ((a >> j) & 1)) {
                hamilton[k][k] += 0.25;
            } else {
                hamilton[k][k] -= 0.25;
                int s = a ^ (1 << i) ^ (1 << j);
                //std::cout << "findstate: ";
                //printBits(s);
                int b = findState(states, s, M);
                hamilton[k][b] = 0.5;
            }
        }
    }
}

// einzelne Blöcke Zusammensetzen
void writeHamiltonBlockToFull(float** hamiltonBlock, float** hamilton, int dimension, int offset) {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            hamilton[i+offset][j+offset] = hamiltonBlock[i][j];
        }
    }
}

// if s is the smallest state integer, returns its periodicity
int checkState(int s, int k) {
    int t = s;
    for (int i = 1; i <= N; i++) {
        t = translateLeft(t, 1);
        if (t < s) {
            return -1;
        } else if (t == s) {
            if (k % (N/i) != 0) {
                return -1;
            } else {
                return i;
            }
        }
    }
    return -1;
}

void checkState(int s, int *r, int *m, int k) {
    //std::cout << "checking state";
    //printBits(s);
    int t = s; *r = -1;
    //std::cout << "def. vars\n";
    for (int i = 1; i <= N; i++) {
        t = translateLeft(t, 1);
        if (t < s) {
            //std::cout << "smaller state found\n";
            return;
        } else if (t == s) {
            //std::cout << "same state found\n";
            if (k % (N/i) != 0) {
                return;
            } else {
                //std::cout << "write to r\n";
                *r = i; continue;
            }
        }
    } //std::cout << "checking state with parity\n";
    t = reflectBits(s); *m = -1;
    for (int i = 0; i < *r; i++) {
        if (t < s) {
            *r = -1; return;
        } else if (t == s) {
            *m = i; return;
        } t = translateLeft(t, 1);
    }
}

void representative(int s, int *r, int *l) {
    int t = s; *r = s; *l = 0;
    for (int i = 1; i < N; i++) {
        //std::cout << "translate\n";
        t = translateLeft(t, 1);
        if (t < *r) {
            *r = t; *l = i;
        }
    }
}

void representative(int s, int *r, int *l, int *q) {
    int t = s; *r = s; *l = 0;
    //std::cout << "vars set";
    for (int i = 1; i < N; i++) {
        //std::cout << "translate\n";
        t = translateLeft(t, 1);
        if (t < *r) {
            *r = t; *l = i;
        }
    } t = reflectBits(s); *q = 0;
    for (int i = 0; i < N; i ++) {
        t = translateLeft(t, 1);
        if (t < *r) {
            *r = t; *l = i; *q = 1;
        }
    }
}

void fillHamiltonMomentumBlock(std::complex<double> **hamilton, std::vector<int> *states, std::vector<int> *R_vals, int k_moment) {
    for (int k = 0; k < states->size(); k++) {
        //std::cout << "a: ";
        int a = states->at(k);
        //printBits(a);
        for (int i = 0; i <= N-1; i++) {
            int j = (i+1) % N;
            if (((a >> i) & 1) == ((a >> j) & 1)) {
                hamilton[k][k] += std::complex<double> (0.25, 0.0);
            } else {
                hamilton[k][k] -= std::complex<double> (0.25, 0.0);
                int s = a ^ (1 << i) ^ (1 << j);
                //printBits(s);
                int r = 0; int l = 0;
                //std::cout << "representative: \n";
                representative(s, &r, &l);
                //std::cout << "findstate: \n";
                int b = findState(states, r, states->size());
                if (b >= 0) {
                    //std::complex<double> numC(0.0, 2 * PI * (double) (bitSum(a) - N/2 ) * (double) l / (double) N );
                    std::complex<double> numC(0.0, 2 * PI * (double) k_moment * (double) l / (double) N );
                    hamilton[k][b] +=  (std::complex<double>) (0.5 * sqrt((double) R_vals->at(k) / (double) R_vals->at(b)) * exp(numC));
                }
            }
        }
    }
}

void writeHamiltonMomentumBlockToFull(std::complex<double> **hamiltonBlock, std::complex<double> **hamilton, int dimension, int offset) {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            //std::cout << i << " " << j << "\n";
            hamilton[i+offset][j+offset] = hamiltonBlock[i][j];
        }
    }
}

void saveComplexHamilton(std::complex<double> **hamilton,const std::string &filename, const std::string &header) {
    std::cout << "saving to file '" << filename << "'..." << std::endl;
    std::ofstream file;
    try {
        file.open(filename);
        file << header << "\n";
        for (int i = 0; i <= size -1; i++) {
            for (int j = 0; j <= size-1; j++) {
                if (hamilton[i][j].real() < 0.001 && hamilton[i][j].real() > -0.001
                && hamilton[i][j].imag() < 0.001 && hamilton[i][j].imag() > -0.001) {
                    file << " \t";
                } else {
                    file << hamilton[i][j] << "\t" ;
                }
            }
            file << std::endl;
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveComplexEiVals(const std::string &filename, const std::string &header, const std::list<std::complex<double>> &eiVals) {
    std::cout << "saving to file '" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open(filename);
        file << header << "\n";
        file << "\nEigenvalues:\n";
        for (std::complex<double> ev : eiVals) {
            file << ev << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void naiverAnatz(double J1, double J2, std::list<std::complex<double>> &H1EiValList) {
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
    saveHamilton(hamilton1, "Hamilton1.txt", "naiver Ansatz für N = " + std::to_string(N));
#endif
#ifdef calculateEigenvalues
    std::cout << "solving...\n";
    Eigen::EigenSolver<Eigen::MatrixXd> solver1(H1);
    const Eigen::VectorXcd &H1EiVal = solver1.eigenvalues();
    for (std::complex<double> ev : H1EiVal) {
        H1EiValList.push_back(ev);
    }
    saveComplexEiVals("EV_naiv.txt", "naiver Ansatz für N = " + std::to_string(N), H1EiValList);
#endif
    for (int i = 0; i < size; i++) {
        delete hamilton1[i];
    } delete[] hamilton1;
}

int main(int argc, char* argv[]) {

    if (argc == 2) {
        std::cout << "N from args\n";
        N = std::stoi(argv[1]);
    }
    std::cout << "no or invalid args given, default N\n";
    size = (int) pow(2, N);
    std::cout << "N: " << N << "; size: " << size << std::endl;

#ifdef naiv
    std::list<std::complex<double>> H1EiValList;
    naiverAnatz(1, 1, H1EiValList);
    std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : H1EiValList) {
        std::cout << ev << "\n";
    }
#endif

    int offset;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef magnetization
    // Using fixed-magnetization blocks.
    std::cout << "\nblockdiagonale m_z" << std::endl;
    static auto **hamilton2 = new float*[size];
    for (int i = 0; i < size; i++) {
        hamilton2[i] = new float[size];
        for (int j = 0; j < size; j++) {
            hamilton2[i][j] = 0.0;
        }
    }

    auto *statesBlock = new std::vector<int>;

    offset = 0;

    for (int n = 0; n <= N; n++) {
        //std::cout << "m_z: " << n << "\n";
        //std::cout << "filling states\n";
        fillStates(statesBlock, n);
        int M = statesBlock->size();
        //std::cout << "statesBlock size: " << M << "\n";
        auto **hamiltonBlock = new float*[M];
        for (int i = 0; i < M; i++) {
            hamiltonBlock[i] = new float[M];
            for (int j = 0; j < M; j++) {
                hamiltonBlock[i][j] = 0.0;
            }
        }
        //std::cout << "filling hamilton\n";
        fillHamiltonBlock(statesBlock, hamiltonBlock, M);
        //std::cout << "writing hamilton block to full\n";
        writeHamiltonBlockToFull(hamiltonBlock, hamilton2, M, offset);
        offset += M;

        statesBlock->clear();
        for (int i = 0; i < M; i++) {
            delete hamiltonBlock[i];
        } delete[] hamiltonBlock;
    }

    Eigen::MatrixXf H2(size, size);
    for (int i = 0; i < size; i++) {
        H2.row(i) = Eigen::VectorXf::Map(&hamilton2[i][0], size);
    }
#ifdef showMatrix
    std::cout << H2 << std::endl;
#endif
#ifdef calculateEigenvalues
    std::cout << "solving...\n";
    Eigen::EigenSolver<Eigen::MatrixXf> solver2(H2);
    Eigen::VectorXcf H2EiVal = solver2.eigenvalues();
    std::cout << H2EiVal << std::endl;
    saveHamilton(hamilton2, "Hamilton2.txt", "Blöcke konstanter Magnetisierung für N = " + std::to_string(N), H2EiVal);
#endif
    for (int i = 0; i < size; i++) {
        delete hamilton2[i];
    } delete[] hamilton2;
#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
