#include <iostream>
#include <fstream>
#include <bitset>
#include <math.h>
#include <string>
#include <vector>
#include <Eigen/Eigenvalues>

#define PI  3.14159265358979323846

// output
#define calculateEigenvalues
#define saveOnlyEiVal
#define showMatrix

// methods
//#define naiv
//#define magnetization
//#define momentum
#define parity

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

// spiegelt die bits
int reflectBits(int s) {
    int t = 0;
    //printBits(s);
    for (int i = 0; i < N; i++) {
        t |= ((s >> (N - 1 - i)) & 1) << i;
    } //printBits(t);
    return t;
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
#ifndef saveOnlyEiVal
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
#endif
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

/////////////////////////???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
double get_gk(int k) {
    if (k == 0 | k == N/2) {
        return 2.0;
    } else {
        return 1.0;
    }
}

/////////////////////////???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
double getNa(int sigma, int m, int Ra, int k, int p) {
    //std::cout << sigma << " " << m << " " << Ra  << " " << k  << " " << p << "\n";
    double Na = pow(N, 2) * get_gk(k) / abs(Ra);
    //std::cout <<  "N: " << Na << "\n";
    if (m != -1) {
        Na *= 1 + sigma * p * cos(k*m);
        //std::cout << "\t" << Na << "\n";
    } return Na;
}

/////////////////////////???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
double helement(int a, int b, int l, int q, int k, int p, std::vector<int> *R_vals, std::vector<int> *m_vals) {
//    if (a == b) {
//        std::cout << "helement same state\n";
//    }
    int sigma_a = R_vals->at(a) / abs(R_vals->at(a));
    int m_a = m_vals->at(a);
    int R_a = R_vals->at(a);
    int sigma_b = R_vals->at(b) / abs(R_vals->at(b));
    int m_b = m_vals->at(b);
    int R_b = R_vals->at(b);
    double Na = getNa(sigma_a, m_a, R_a, k, p);
    double Nb = getNa(sigma_b, m_b, R_b, k, p);

    std::cout << "\t" << k << " " << l << " " << sigma_a << " " << sigma_b<< " " << p << " " << m_b << "\n";

    double val = 0.5 * (double) pow(sigma_a * p, q) * sqrt(Na / Nb);// <--------- für h_j(a) = 0.5 einsetzen????
    //std::cout << val << "\n";
    //std::cout << "\t" << Na << " " << Nb << " " << sigma_a << " " << p << " " << q << " " << pow(sigma_a * p, q) << "\n";

    if (sigma_a == sigma_b) { // same sigma
        if (m_b != -1) {
            std::cout << "same sigma, m != -1, changed from " << val;
            val *= (cos(k*l) + sigma_a * p * cos(k*(l-m_b))) / (1 + sigma_a * p * cos(k*m_b));
            std::cout << " to " << val << " by multiplaying ";
            std::cout << (cos(k*l) + sigma_a * p * cos(k*(l-m_b))) / (1 + sigma_a * p * cos(k*m_b));
        } else {
            std::cout << "same sigma changed from " << val;
            val *= cos(k*l);
            std::cout << " to " << val << " by multiplaying " << cos(k*l) << "\n";
        }
    } else {
        if (m_b != -1) {
            std::cout << "different sigma, m != -1, changed from " << val;
            val *= (-sigma_a * sin(k*l) + p * sin(k*(l-m_b))) / (1 - sigma_a * p * cos(k*m_b));
            std::cout << " to " << val << " by multiplaying " << (-sigma_a * sin(k*l) + p * sin(k*(l-m_b))) / (1 - sigma_a * p * cos(k*m_b)) << "\n";
        } else {
            std::cout << "different sigma changed from " << val;
            val *= -sigma_a * sin(k*l);
            std::cout << " to " << val << " by multiplaying " << -sigma_a * sin(k*l) << "\n";
        }
    }
    std::cout << "helement returned: " << val << "\n";
    return val;
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
#ifndef saveOnlyEiVal
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
#endif
        file << "\nEigenvalues:\n";
        file << eiVal;
    } catch (...) {
        std::cout << "failed to save to file\n";
    }

    std::cout << "done\n";
    file.close();

}

void fillHamiltonParityBlock(std::complex<double> **hamilton, std::vector<int> *states, std::vector<int> *R_vals, std::vector<int> *m_vals, int k_moment, int p) {
    for (int k = 0; k < states->size(); k++) {
        //std::cout << "setting n\n";
        int n = 1;
        if (k > 0 && states->at(k-1) == states->at(k)) {
            continue;
        } else if (k < states->size() - 1 && states->at(k) == states->at(k+1)) {
            n = 2;
        }

        //std::cout << "getting state\n";
        int a = states->at(k);
        std::complex<double> Ez(0.0, 0.0);

        for (int i = 0; i <= N-1; i++) {
            //std::cout << "loop" << i << "\n";
            int j = (i+1) % N;
            if (((a >> i) & 1) == ((a >> j) & 1)) {
                Ez += std::complex<double> (0.25, 0.0);
            } else {
                Ez -= std::complex<double> (0.25, 0.0);
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                int s = a ^ (1 << i) ^ (1 << j);
                int r = 0, l = 0, q = 0;
                //std::cout << "representative\n";
                representative(s, &r, &l, &q);
                //std::cout << "findState\n";
                int b = findState(states, r, states->size());
                std::cout << k_moment << " " << r << " " << b << "\n";
                int m;
                if (b >= 0) {
                    //std::cout << "setting m\n";
                    if (b > 1 && states->at(b) == states->at(b - 1)) {
                        m = 2, b = b - 1;
                    } else if (b < states->size() - 1 && states->at(b) == states->at(b + 1)) {
                        m = 2;
                    } else {
                        m = 1;
                    }
                    for (int nested_j = b; nested_j < b + m; nested_j++) {
                        for (int nested_i = k; nested_i < k + n; nested_i++) {
                            //std::cout << "changed element " << nested_i << " " << nested_j << "from " << hamilton[nested_i][nested_j];
                            hamilton[nested_i][nested_j] += (std::complex<double>) helement(nested_i, nested_j, l, q, k_moment, p, R_vals, m_vals);
                            //std::cout << " to " << hamilton[nested_i][nested_j] << "\n";
                            //std::cout << nested_i << " " << nested_j << " " << (std::complex<double>) helement(nested_i, nested_j, l, q, k, p, R_vals, m_vals) << "\n";
                        }
                    }
                }
            }
        }
        //std::cout << hamilton[k][k] << " " << n << " " << Ez << "\n";
        hamilton[k][k] += (double) n * Ez;
    }
}

int main(int argc, char* argv[]) {

    std::cout << "N: " << N << "; size: " << size << std::endl;

#ifdef naiv
    // Methode 1
    std::cout << "\nnaiver Ansatz:..." << std::endl;
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
#ifdef showMatrix
    std::cout << H1 << std::endl;
#endif
#ifdef calculateEigenvalues
    std::cout << "solving...\n";
    Eigen::EigenSolver<Eigen::MatrixXf> solver1(H1);
    Eigen::VectorXcf H1EiVal = solver1.eigenvalues();
    std::cout << H1EiVal << std::endl;
    saveHamilton(hamilton1, "Hamilton1.txt", "naiver Ansatz für N = " + std::to_string(N), H1EiVal);
#endif
    for (int i = 0; i < size; i++) {
        delete hamilton1[i];
    } delete[] hamilton1;
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
        std::cout << "k: " << k << "\n";
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
                    printBits(a);
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

#ifdef  parity
    // Using the parity operator reflectBits()
    std::cout << "\nparity operator:..." << std::endl;

    //std::cout << "allocating matrix\n";
    static auto **hamilton4 = new std::complex<double>*[size];
    for (int i = 0; i < size; i++) {
        hamilton4[i] = new std::complex<double>[size];
        for (int j = 0; j < size; j++) {
            hamilton4[i][j] = std::complex<double> (0.0, 0.0);
        }
    }


    auto *statesP = new std::vector<int>;
    auto *statesListP = new std::vector<int>;
    auto *statesPerioP = new std::vector<int>;
    auto *statesTransP = new std::vector<int>;
    offset = 0;
    for (int k = -(N-1)/2; k <= N/2 ; k++) { // go though momentums
        std::cout << "momentum k: " << k << "\n";
        for (int m = 0; m <= N; m++) { // go though magnetizations
            //std::cout << "magnetization m: " << m << "\n";
            fillStates(statesP, m);
            for (int a : *statesP) {
                //std::cout << "found state: ";
                //printBits(a);
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                int perio = 0, m_trans = 0;
                checkState(a, &perio, &m_trans, k);
                //std::cout << perio << " " << m_trans << " \n";
                //printBits(a);
                //std::cout << "checked state\n";
                for (int p : {-1, 1}) {
                    if ((k != 0) && (k != N/2) && p == -1) {
                        continue;
                    }
                    for (int sigma : {-1, 1}) {
                        if ((k != 0) && (k != N/2) && sigma == 1) {
                            continue;
                        }
                        if (m_trans != -1) {
                            std::complex<double> cosArg(0.0, k * m_trans * 2 * PI / N);
                            std::complex<double> val = (double) sigma * (double) p * cos(cosArg);
                            //std::cout << sigma << " " << cosArg << " " << val << "\n";
                            //std::cout << "val: " << val << "\n";
                            if ( (abs((double) 1 + val) < 0.001) | (sigma == -1 && abs((double) 1 - val) > 0.001)) {
                                //std::cout << "invalid\n";
                                perio = -1;
                            }
                        } if (perio > 0) {
                            //std::cout << "parity p: " << p << "\n";
                            //std::cout << "sigma: " << sigma << "\n";
                            //std::cout << "m_trans: " << m_trans << "\n";
                            std::cout << "valid state: ";
                            printBits(a);
                            //std::cout << "adding to statesList\n";
                            statesListP->push_back(a);
                            statesPerioP->push_back(sigma * perio);
                            statesTransP->push_back(m_trans);
                        }
                    }
                    int M = statesListP->size();
                    auto **hamiltonBlock = new std::complex<double>*[M];
                    for (int i = 0; i < M; i++) {
                        hamiltonBlock[i] = new std::complex<double>[M];
                        for (int j = 0; j < M; j++) {
                            hamiltonBlock[i][j] = std::complex<double> (0.0, 0.0);
                        }
                    }
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //std::cout << "filling hamilton block" << M << "\n";
                    fillHamiltonParityBlock(hamiltonBlock, statesListP, statesPerioP, statesTransP, k, p);
                    //std::cout << "writing hamilton block to full\n";
                    writeHamiltonMomentumBlockToFull(hamiltonBlock, hamilton4, M, offset);
                    offset += M;
                    //std::cout << offset << "\n";
                    //std::cout << "clean up\n";
                    statesP->clear();
                    statesListP->clear();
                    statesPerioP->clear();
                    statesTransP->clear();
                    for (int i = 0; i < M; i++) {
                        delete hamiltonBlock[i];
                    } delete[] hamiltonBlock;
                }
            }
        }
    }

    Eigen::MatrixXcd H4(size, size);
    for (int i = 0; i < size; i++) {
        H4.row(i) = Eigen::RowVectorXcd::Map(&hamilton4[i][0], size);
    }
#ifdef showMatrix
    std::cout << H4 << std::endl;
#endif
#ifdef calculateEigenvalues
    std::cout << "solving...\n";
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver4(H4);
    Eigen::VectorXcd H4EiVal = solver4.eigenvalues();
    std::cout << H4EiVal << std::endl;
    saveComplexHamilton(hamilton4, "Hamilton4.txt", "momentum states", H4EiVal);
#endif
    for (int i = 0; i < size; i++) {
        delete hamilton4[i];
    } delete[] hamilton4;
#endif

    return 0;
}
