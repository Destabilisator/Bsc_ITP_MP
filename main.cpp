#include <iostream>
#include <fstream>
#include <bitset>
#include <math.h>
#include <string>
#include <vector>
#include <Eigen/Eigenvalues>

#define naiv
//#define magnetization
//#define momentum

const short N = 8;
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

void saveHamilton(float** hamilton, std::string filename, std::string header) {
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
                    file <<  hamilton[i][j] << "\t" ;
                }
            }
            file << std::endl;
        }

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

// findet alle Zustände mit der Magnetisierung m_z = n
void fillStates(std::vector<int> *states, int n) {
    for (int s = 0; s <= size - 1; s++) {
        if (bitSum(s) == n) {
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
                int b = a ^ (1 << i) ^ (1 << j);
                //std::cout << "findstate: ";
                //printBits(b);
                int pos = findState(states, b, M);
                hamilton[k][pos] = 0.5;
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

int main(int argc, char* argv[]) {

    std::cout << "N: " << N << "; size: " << size << std::endl;

    // Methode 1
#ifdef naiv
    std::cout << "naiver Ansatz" << std::endl;
    static auto **hamilton1 = new float*[size];
    for (int i = 0; i < size; i++) {
        hamilton1[i] = new float[size];
        for (int j = 0; j < size; j++) {
            hamilton1[i][j] = 0.0;
        }
    }

    fillHamilton1(hamilton1);
    //saveHamilton(hamilton1, "Hamilton1.txt", "naiver Ansatz für N = " + std::to_string(N));

    Eigen::MatrixXf H1(size, size);
    for (int i = 0; i < size; i++) {
        H1.row(i) = Eigen::VectorXf::Map(&hamilton1[i][0], size);
    }
    //std::cout << H1 << std::endl;
    std::cout << "solving... ";
    Eigen::EigenSolver<Eigen::MatrixXf> solver1(H1);
    std::cout << "done\n";
    std::cout << solver1.eigenvalues() << std::endl;

#endif

#ifdef magnetization
    // Using fixed-magnetization blocks.
    std::cout << "blockdiagonale m_z" << std::endl;

    static auto **hamilton2 = new float*[size];
    for (int i = 0; i < size; i++) {
        hamilton2[i] = new float[size];
        for (int j = 0; j < size; j++) {
            hamilton2[i][j] = 0.0;
        }
    }

    auto *statesBlock = new std::vector<int>;

    int offset = 0;

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
        writeHamiltonBlockToFull(hamiltonBlock, hamilton2, M, offset);
        offset += M;

        statesBlock->clear();
        for (int i = 0; i < M; i++) {
            delete hamiltonBlock[i];
        } delete[] hamiltonBlock;
    }

    saveHamilton(hamilton2, "Hamilton2.txt", "Blöcke konstanter Magnetisierung für N = " + std::to_string(N));

    Eigen::MatrixXf H2(size, size);
    for (int i = 0; i < size; i++) {
        H2.row(i) = Eigen::VectorXf::Map(&hamilton2[i][0], size);
    }
    std::cout << H2 << std::endl;
    Eigen::EigenSolver<Eigen::MatrixXf> solver2(H2);
    Eigen::VectorXf H2EiVal = solver2.eigenvalues();
    std::cout << H2EiVal << std::endl;


#endif

#ifdef momentum
    auto *statesList = new std::vector<int>;
    auto *states = new std::vector<int>;
    auto *statesPerio = new std::vector<int>;
    for (int k = 0; k <= N; k++) {
        fillStates(states, k);
        for (int i = 0; i < states->size(); i++) {
            int a = states->at(i);
            //std::cout << "found state: ";
            //printBits(a);
            int perio = checkState(a, k);
            if (perio >= 0) {
                //std::cout << perio << "\n";
                statesList->push_back(a);
                //std::cout << statesList->at(0) << ": ";
                statesPerio->push_back(perio);
                //std::cout << statesPerio->at(0) << "\n";
            }
        } states->clear();
    }

    for (int i = 0; i < statesList->size(); i++) {
        std::cout << statesPerio->at(i) << ": ";
        printBits(statesList->at(i));
    }

#endif

    return 0;
}
