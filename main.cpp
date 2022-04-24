#include <iostream>
#include <fstream>
#include <bitset>
#include <math.h>
#include <string>
#include <vector>

const short N = 3;
int size = (int) pow(2, N);

void printBits(int a) {
    std::bitset<N> x(a);
    std::cout << x << '\n';
}

int bitSum(int s) {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        sum += ((s >> i) & 1 );
    } return sum;
}

void saveHamilton(float** hamilton, std::string filename) {
    std::cout << "saving to file '" << filename << "'..." << std::endl;

    std::ofstream file;
    try {
        file.open(filename);
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

void fillHamilton1(float** hamilton) {
    for (int a = 0; a <= size -1; a++) {
        std::cout << "a: ";
        printBits(a);
        for (int i = 0; i <= N-1; i++) {
            int j = (i+1) % N;
            if (((a >> i) & 1) == ((a >> j) & 1)) {
                std::cout << "on " <<  a << "\n";
                hamilton[a][a] += 0.25;
            } else {
                hamilton[a][a] -= 0.25;
                int b = a ^ (1 << i) ^ (1 << j);
                std::cout << "off " <<  a << " " << b << "\n";
                hamilton[a][b] = 0.5;
            }
        }
    }
}

void fillStates(std::vector<int> *states, int n) {
    for (int s = 0; s <= size - 1; s++) {
        if (bitSum(s) == n) {
            //std::cout << bitSum(s) << ": ";
            states->push_back(s);
            //printBits(s);
        }
    } states->shrink_to_fit();
}


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

void writeHamiltonBlockToFull(float** hamiltonBlock, float** hamilton, int dimension, int offset) {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            hamilton[i+offset][j+offset] = hamiltonBlock[i][j];
        }
    }
}

int main() {
    std::cout << "N: " << N << "; size: " << size << std::endl;

    // Methode 1

    std::cout << "declaring array" << std::endl;
    static auto **hamilton1 = new float*[size];
    for (int i = 0; i < size; i++) {
        hamilton1[i] = new float[size];
        for (int j = 0; j < size; j++) {
            hamilton1[i][j] = 0.0;
        }
    }

    std::cout << "writing to array" << std::endl;
    fillHamilton1(hamilton1);
    saveHamilton(hamilton1, "Hamilton1.txt");

    // Using fixed-magnetization blocks.

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

    saveHamilton(hamilton2, "Hamilton2.txt");

    return 0;
}
