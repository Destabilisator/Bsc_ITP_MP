#include <iostream>
#include <fstream>
#include <bitset>
#include <math.h>
#include <string>
#include <vector>
#include <Eigen/Eigenvalues>

#define PI  3.14159265358979323846

#define calculateEigenvalues

#define naiv
#define magnetization
#define momentum

const short N = 3;
const int size = (int) pow(2, N);

void printBits(int a) {
    std::bitset<N> x(a);
    std::cout << x << '\n';
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

// sum up all bits in s
int bitSum(int s) {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        sum += ((s >> i) & 1 );
    } return sum;
}

void saveHamilton(float** hamilton, std::string filename, std::string header, Eigen::VectorXcf eiVal) {
    std::cout << "saving to file '" << filename << "'..." << std::endl;

    std::ofstream file;
    try {
        file.open(filename);
        file << header << "\n";
        for (int i = 0; i <= size -1; i++) {
            for (int j = 0; j <= size-1; j++) {
                if (hamilton[i][j] < 0.1 && hamilton[i][j] > -0.1) {
                    file << " \t";
                } else {
                    file << hamilton[i][j] << "\t" ;
                }
            }
            file << std::endl;
        }
        file << "\nEigenvalues:\n";
        file << eiVal;
    } catch (...) {
        std::cout << "failed to save to file\n";
    }

    std::cout << "done\n";
    file.close();

}

// naiver Ansatz: direkt füllen, ohne zu ordnen
void fillHamilton1(float** hamilton) {
    for (int a = 0; a <= size -1; a++) {
        //std::cout << "a: ";
        //printBits(a);
        for (int i = 0; i <= N-1; i++) {
            int j = (i+1) % N;
            if (((a >> i) & 1) == ((a >> j) & 1)) {
                //std::cout << "on " <<  a << "\n";
                hamilton[a][a] += 0.25;
            } else {
                hamilton[a][a] -= 0.25;
                int b = a ^ (1 << i) ^ (1 << j);
                //std::cout << "off " <<  a << " " << b << "\n";
                hamilton[a][b] = 0.5;
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
    for (int i = 1; i <= N + 1; i++) {
        t = translateLeft(t, 1);
        //printBits(t);
        if (t < s) {
            //std::cout << "found smaller\n";
            return -1;
        } else if (t == s) {
            if (k % (N/i) != 0) {
                //std::cout << "not compatible\n";
                return -1;
            } else {
                //std::cout << "found compatible state\n";
                return i;
            }
        }
    }
    return -1;
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

void saveComplexHamilton(std::complex<double> **hamilton, std::string filename, std::string header, Eigen::VectorXcd eiVal) {
    std::cout << "saving to file '" << filename << "'..." << std::endl;

    std::ofstream file;
    try {
        file.open(filename);
        file << header << "\n";
        for (int i = 0; i <= size -1; i++) {
            for (int j = 0; j <= size-1; j++) {
                if (hamilton[i][j].real() < 0.1 && hamilton[i][j].real() > -0.1
                && hamilton[i][j].imag() < 0.1 && hamilton[i][j].imag() > -0.1) {
                    file << " \t";
                } else {
                    file << hamilton[i][j] << "\t" ;
                }
            }
            file << std::endl;
        }
        file << "\nEigenvalues:\n";
        file << eiVal;
    } catch (...) {
        std::cout << "failed to save to file\n";
    }

    std::cout << "done\n";
    file.close();

}

int main(int argc, char* argv[]) {

    std::cout << "N: " << N << "; size: " << size << std::endl;

#ifdef naiv
    // Methode 1
    std::cout << "\nnaiver Ansatz" << std::endl;
    static auto **hamilton1 = new float*[size];
    for (int i = 0; i < size; i++) {
        hamilton1[i] = new float[size];
        for (int j = 0; j < size; j++) {
            hamilton1[i][j] = 0.0;
        }
    }

    fillHamilton1(hamilton1);

    Eigen::MatrixXf H1(size, size);
    for (int i = 0; i < size; i++) {
        H1.row(i) = Eigen::VectorXf::Map(&hamilton1[i][0], size);
    }
    std::cout << H1 << std::endl;
#ifdef calculateEigenvalues
    std::cout << "solving...";
    Eigen::EigenSolver<Eigen::MatrixXf> solver1(H1);
    Eigen::VectorXcf H1EiVal = solver1.eigenvalues();
    std::cout << H1EiVal << std::endl;
    saveHamilton(hamilton1, "Hamilton1.txt", "naiver Ansatz für N = " + std::to_string(N), H1EiVal);

#endif

#endif

    int offset;

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
    std::cout << H2 << std::endl;
#ifdef calculateEigenvalues
    std::cout << "solving...\n";
    Eigen::EigenSolver<Eigen::MatrixXf> solver2(H2);
    Eigen::VectorXcf H2EiVal = solver2.eigenvalues();
    std::cout << H2EiVal << std::endl;
    saveHamilton(hamilton2, "Hamilton2.txt", "Blöcke konstanter Magnetisierung für N = " + std::to_string(N), H2EiVal);
#endif

#endif

#ifdef momentum
    // Using momentum states.
    std::cout << "\nmomentum states" << std::endl;

    //std::cout << "allocating matrix\n";
    static auto **hamilton3 = new std::complex<double>*[size];
    for (int i = 0; i < size; i++) {
        hamilton3[i] = new std::complex<double>[size];
        for (int j = 0; j < size; j++) {
            hamilton3[i][j] = std::complex<double> (0.0, 0.0);
        }
    }


    auto *states = new std::vector<int>;
    auto *statesList = new std::vector<int>;
    auto *statesPerio = new std::vector<int>;
    offset = 0;
    for (int k = -N/2; k <= N/2 ; k++) { // <---------------------------------------------------------------- Problemski
        //std::cout << "k: " << k << "\n";
        for (int m = 0; m <= N; m++) {
            fillStates(states, m);
            for (int a : *states) {
                //std::cout << "found state: ";
                //printBits(a);
                int perio = checkState(a, k);
                if (perio >= 0) {
                    //std::cout << perio << "\n";
                    statesList->push_back(a);
                    //std::cout << statesList->at(0) << ": ";
                    statesPerio->push_back(perio);
                    //printBits(a);
                    //std::cout << statesPerio->at(0) << "\n";
                }
            }
            int M = statesList->size();
            auto **hamiltonBlock = new std::complex<double>*[M];
            for (int i = 0; i < M; i++) {
                hamiltonBlock[i] = new std::complex<double>[M];
                for (int j = 0; j < M; j++) {
                    hamiltonBlock[i][j] = std::complex<double> (0.0, 0.0);
                }
            }
            //std::cout << "filling hamilton block\n";
            fillHamiltonMomentumBlock(hamiltonBlock, statesList, statesPerio, k);
            //std::cout << "writing hamilton block to full\n";
            writeHamiltonMomentumBlockToFull(hamiltonBlock, hamilton3, M, offset);
            offset += M;

            //std::cout << "clean up\n";
            states->clear();
            statesList->clear();
            statesPerio->clear();
            for (int i = 0; i < M; i++) {
                delete hamiltonBlock[i];
            } delete[] hamiltonBlock;
        }

    }

    Eigen::MatrixXcd H3(size, size);
    for (int i = 0; i < size; i++) {
        H3.row(i) = Eigen::RowVectorXcd::Map(&hamilton3[i][0], size);
    }
    std::cout << H3 << std::endl;
#ifdef calculateEigenvalues
    std::cout << "solving...\n";
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver3(H3);
    Eigen::VectorXcd H3EiVal = solver3.eigenvalues();
    std::cout << H3EiVal << std::endl;
    saveComplexHamilton(hamilton3, "Hamilton3.txt", "momentum states", H3EiVal);
#endif

#endif

    return 0;
}
