#include "helpers.h"

/////////////////////////////// bits ///////////////////////////////

void printBits(int s, int N) {
    for (int n = N-1; n >= 0; n--) {
        std::cout << ( (s >> n) & 1 );
    } std::cout << "\n";
}

int translateRight(int s, int n, int N) {
    for (int _ = 0; _ < n; _++) {
        int bit = s & 1;
        s = (s>>1) | (bit << (N-1));
    } return s;
}

int translateLeft(int s, int n, int N) {
    for (int _ = 0; _ < n; _++) {
        int bit = (s >> (N-1)) & 1;
        s = (s<<1) | bit;
        s &= ~(1 << N);
    } return s;
}

int reflectBits(int s, int N) {
    int t = 0;
    for (int i = 0; i < N; i++) {
        t |= ((s >> (N - 1 - i)) & 1) << i;
    }
    return t;
}

int invertBits(int s, int N) {
    return s ^ (int) pow(2,N)-1;
}

int bitSum(int s, int N) {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        sum += ((s >> i) & 1 );
    } return sum;
}

/////////////////////////////// states ///////////////////////////////

void fillStates(std::vector<int> *states, int m, int N, int size) {
    for (int s = 0; s <= size - 1; s++) {
        if (bitSum(s, N) == m) {
            states->push_back(s);
        }
    } states->shrink_to_fit();
}

int findState(const std::vector<int>& states, int s) {
    int pos, pos_min = 0, pos_max = states.size()-1;
    while (true) {
        pos = pos_min + (pos_max - pos_min ) / 2;
        if (s < states.at(pos)) {
            pos_max = pos - 1;
        } else if (s > states.at(pos)) {
            pos_min = pos + 1;
        } else {
            return pos;
        } if (pos_min > pos_max) {
            return -1;
        }
    }
}

int checkState(int s, int k, int N) {
    int t = s;
    for (int i = 1; i <= N; i++) {
        t = translateLeft(t, 1, N);
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

void checkState(int s, int *r, int *m, int k, int N) {
    int t = s; *r = -1;
    for (int i = 1; i <= N; i++) {
        t = translateLeft(t, 1, N);
        if (t < s) {
            return;
        } else if (t == s) {
            if (k % (N / i) != 0) {
                return;
            } else {
                *r = i;
                continue;
            }
        }
    }
    t = reflectBits(s, N); *m = -1;
    for (int i = 0; i < *r; i++) {
        if (t < s) {
            *r = -1; return;
        } else if (t == s) {
            *m = i; return;
        } t = translateLeft(t, 1, N);
    }
}

void representative(int s, int *r, int *l, int N) {
    int t = s; *r = s; *l = 0;
    for (int i = 1; i < N; i++) {
        t = translateLeft(t, 1, N);
        if (t < *r) {
            *r = t; *l = i;
        }
    }
}

void representative(int s, int *r, int *l, int *q, int N) {
    int t = s; *r = s; *l = 0;
    for (int i = 1; i < N; i++) {
        t = translateLeft(t, 1, N);
        if (t < *r) {
            *r = t; *l = i;
        }
    } t = reflectBits(s, N); *q = 0;
    for (int i = 0; i < N; i ++) {
        t = translateLeft(t, 1, N);
        if (t < *r) {
            *r = t; *l = i; *q = 1;
        }
    }
}

/////////////////////////////// saving data ///////////////////////////////

void saveEiVals(const std::string &filename, const std::string &header, const std::list<double> &eiVals) {
    std::cout << "saving to file '" << filename << "'..." << std::endl;
    std::ofstream file;
    try {
        file.open(filename);
        file << header <<".\n";
        file << "Eigenvalues:\n";
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

void saveComplexEiVals(const std::string &filename, const std::string &header, const std::list<std::complex<double>> &eiVals) {
    std::cout << "saving to file '" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open(filename);
        file << header << ".\n";
        file << "Eigenvalues:\n";
        for (std::complex<double> ev : eiVals) {
            file << ev << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveHamilton(double** hamilton, const std::string &filename, const std::string &header, int size) {
    std::cout << "saving to file '" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open(filename);
        file << header << "\n";
        for (int i = 0; i <= size -1; i++) {
            for (int j = 0; j <= size-1; j++) {
                //if (hamilton[i][j] < 0.001 && hamilton[i][j] > -0.001) {
                //    file << " \t";
                //} else {
                    file << hamilton[i][j] << "\t" ;
                //}
            }
            file << "\n";
        }

    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveComplexHamilton(std::complex<double> **hamilton,const std::string &filename, const std::string &header, int size) {
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
            file << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}
