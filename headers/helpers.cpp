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
    int pos, pos_min = 0, pos_max = (int) states.size()-1;
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

int findState(const std::vector<std::tuple<int, int>>& states, int s) {
    int pos, pos_min = 0, pos_max = (int) states.size()-1;
    while (true) {
        pos = pos_min + (pos_max - pos_min ) / 2;
        if (s < std::get<0>(states.at(pos))) {
            pos_max = pos - 1;
        } else if (s > std::get<0>(states.at(pos))) {
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
    for (int i = 1; i <= N/2; i++) {
        t = translateLeft(t, 2, N);
        if (t < s) {return -1;} else if (t == s) {
            if (k % (int) ((double) N / (double) i / 2.0) != 0) {return -1;}
            else {return i;}
        }
    }
    return -1;
}

void checkState(int s, int *r, int *m, int k, int N) {
    int t = s; *r = -1; *m = -1;
    //for (int i = 1; i <= N; i++) {
    for (int i = 1; i <= N/2; i++) {
        //t = translateLeft(t, 1, N);
        t = translateLeft(t, 2, N);
        if (t < s) {return;}
        else if (t == s) {
            //if (k % (int) ((double) N / (double) i) != 0) {return;}
            if (k % (int) ((double) N / (double) i / 2.0) != 0) {return;}
            else {*r = i; break;}
        }
    }
    //t = reflectBits(s, N);
    t = reflectBits(translateLeft(s, 1, N), N);
    //t = translateRight(t, 1, N);
    for (int i = 0; i < *r; i++) {
        if (t < s) {*r = -1; return;}
        else if (t == s) {*m = i; return;}
        //t = translateLeft(t, 1, N);
        t = translateLeft(t, 2, N);
    }
}

void representative(int s, int *r, int *l, int N) {
    int t = s; *r = s; *l = 0;
    for (int i = 1; i < N/2; i++) {
        t = translateLeft(t, 2, N);
        if (t < *r) {*r = t; *l = i;}
    }
}

void representative(int s, int *r, int *l, int *q, int N) {
    int t = s; *r = s; *l = 0;
    for (int i = 1; i < N/2; i++) {
        t = translateLeft(t, 2, N);
        if (t < *r) {*r = t; *l = i;}
    }
    t = reflectBits(s, N); *q = 0;
    for (int i = 1; i < N/2; i ++) {
        t = translateLeft(t, 2, N);
        if (t < *r) {*r = t; *l = i; *q = 1;}
    }
}

/////////////////////////////// saving data ///////////////////////////////

void saveEiVals(const std::string &filename, const std::string &header, const std::list<double> &eiVals, const int &N) {
    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << std::endl;
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header <<"\n\n";
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

void saveEiVals(const std::string &filename, const std::string &header, const std::vector<double> &eiVals, const int &N) {
    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
        file << "Eigenvalues:\n";
        for (double ev : eiVals) {
            file << ev << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveComplexEiVals(const std::string &filename, const std::string &header, const std::list<std::complex<double>> &eiVals, const int &N) {
    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
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

void saveComplexEiVals(const std::string &filename, const std::string &header, const std::vector<std::complex<double>> &eiVals, const int &N) {
    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
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

void saveHamilton(double** hamilton, const std::string &filename, const std::string &header, const int &size, const int &N) {
    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
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

void saveComplexHamilton(std::complex<double> **hamilton,const std::string &filename, const std::string &header, const int &size, const int &N) {
    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << std::endl;
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
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

void saveMatrixToFile(const Eigen::MatrixXd& matrix, const std::string &filename, const std::string &header, const int &N) {
    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
        file << matrix;
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveComplexMatrixToFile(const Eigen::MatrixXcd& matrix, const std::string &filename, const std::string &header, const int &N) {
    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
        file << matrix;
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveOutData(const std::string &filename, const std::string &header, const std::string &x_label,
                 const std::string &y_label, const std::vector<std::tuple<double, double>> &outData, const int &N) {

    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
        file << x_label << "\t" << y_label << "\n";
        for (std::tuple<double, double> data : outData) {
            file << std::get<0>(data) << "\t" << std::get<1>(data) << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveOutData(const std::string &filename, const std::string &header, const std::string &x_label,
                 const std::string &y_label, const std::vector<std::tuple<int, double>> &outData, const int &N) {

    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
        file << x_label << "\t" << y_label << "\n";
        for (std::tuple<double, double> data : outData) {
            file << std::get<0>(data) << "\t" << std::get<1>(data) << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

/////////////////////////////// calculate quantities ///////////////////////////////

Eigen::MatrixXd spinMatrix(const int &N, const  int &SIZE) {
    Eigen::MatrixXd S2 = 0.75 * (double) N * Eigen::MatrixXd::Identity(SIZE, SIZE);
    for (int s = 0; s < SIZE; s++) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < j; i++) {
                if (((s >> i) & 1) == ((s >> j) & 1)) {
                    S2(s, s) += 0.5;
                } else {
                    S2(s, s) -= 0.5;
                    int d = s ^ (1 << i) ^ (1 << j);
                    S2(s, d) = 1.0;
                }
            }
        }
    }
    return S2;
}

Eigen::MatrixXd spinMatrix(const int &N, const std::vector<int> &states) {
    int size = (int) states.size();
    Eigen::MatrixXd S2 = 0.75 * (double) N * Eigen::MatrixXd::Identity(size, size);
    for (int k = 0; k < size; k++) {
        for (int j = 0; j < N; j++) {
            int s = states.at(k);
            for (int i = 0; i < j; i++) {
                if (((s >> i) & 1) == ((s >> j) & 1)) {
                    S2(k, k) += 0.5;
                } else {
                    S2(k, k) -= 0.5;
                    int d = s ^ (1 << i) ^ (1 << j);
                    int b = findState(states, d);
                    S2(k, b) = 1.0;
                }
            }
        }
    }
    return S2;
}

double getSpecificHeat(const double &temp, const std::vector<std::complex<double>>& eiVals, const int &N) {
    double Z_sum = 0.0, expectation_H = 0.0, expectation_H_2 = 0.0;
    for (std::complex<double> ev : eiVals) {
        double ev_real = std::real(ev);
        Z_sum += std::exp(-1.0 / temp * ev_real);
        expectation_H += std::exp(-1.0 / temp * ev_real) * ev_real;
        expectation_H_2 += std::exp(-1.0 / temp * ev_real) * ev_real * ev_real;
    }
    expectation_H /= Z_sum;
    expectation_H_2 /= Z_sum;
    return 1.0 / temp * 1.0 / temp * ( expectation_H_2 - expectation_H * expectation_H ) / N;
}

double getSusceptibility(const double &temp, const Eigen::MatrixXcd &M, const std::vector<std::complex<double>>& eiVals, const int &N) {
    double Z_sum = 0.0, expectation_mz_2 = 0.0;
    for (int i = 0; i < eiVals.size(); i++) {
        double ev_real = std::real(eiVals.at(i));
        Z_sum += std::exp(-1.0 / temp * ev_real);
        expectation_mz_2 += std::exp(-1.0 / temp * ev_real) * std::real(M(i, i));
    }
    expectation_mz_2 /= Z_sum;
    expectation_mz_2 /= 3.0;
    return 1.0 / temp * expectation_mz_2 / N;
}

double getSusceptibilityDegeneracy(const double &temp, const Eigen::MatrixXcd &M, const std::vector<std::complex<double>>& eiVals, const int &N) {
    double Z_sum = 0.0, expectation_mz_2 = 0.0;
    for (int i = 0; i < eiVals.size(); i++) {
        double ev_real = std::real(eiVals.at(i));
        double S_elem = std::real(M(i, i));
        double S = - 0.5 + std::sqrt(0.25 + S_elem);
        Z_sum += std::exp(-1.0 / temp * ev_real) * (2.0 * S + 1);
        expectation_mz_2 += std::exp(-1.0 / temp * ev_real) * S_elem * (2.0 * S + 1);
    }
    expectation_mz_2 /= Z_sum;
    expectation_mz_2 /= 3.0;
    return 1.0 / temp * expectation_mz_2 / N;
}

/////////////////////////////// others ///////////////////////////////

// [executable] N J_START J_END J_COUNT CORES SILENT
void validateInput(int argc, char* argv[], int *N, int *SIZE, double *J_START, double *J_END, int *J_COUNT,
                   const unsigned int *cpu_cnt, bool *silent, int *cores, const double *J1, const double *J2) {

    if (argc >= 2) {
        if ( (std::stoi(argv[1]) % 2 == 0) && (std::stoi(argv[1]) >= 6) && (std::stoi(argv[1]) <= 32) ) {
            *N = std::stoi(argv[1]);
        } else {
            std::cout << "invalid chain size, must be even and at least 6, defaulting to " << N << "\n";
        }
    }
    *SIZE = (int) pow(2, *N);
    std::cout << "N: " << *N << "; size: " << *SIZE << std::endl;
    if (argc >= 5) {
        std::cout << "range given: ";
        if (std::stod(argv[2]) > std::stod(argv[3]) || std::stoi(argv[4]) < 1) {
            std::cout << "range invalid, defaulting...\n";
        } else {
            *J_START = std::stod(argv[2]); *J_END = std::stod(argv[3]); *J_COUNT = std::stoi(argv[4]);
        }
        std::cout << "START = " << *J_START << ", END = " << *J_END << " and COUNT = " << *J_COUNT << " from args\n";
    } else {
        std::cout << "no range given: ";
        std::cout << "START = " << *J_START << ", END = " << *J_END << " and COUNT = " << *J_COUNT << " from default\n";
    }

    std::cout << "default J1 and J2 for beta plot: J1 = " << *J1 << " and J2 = " << *J2 << " (currently unchangeable)\n";

    if (argc >= 6) {
        int crs = std::stoi(argv[5]);
        if (crs > 0 && crs <= *cpu_cnt) {
            *cores = crs;
            std::cout << "using " << *cores << " cores\n";
        } else {
            std::cout << "defaulting to using all (" << *cores << ") cores\n";
        }
    }

    return;

    if (argc >= 7) {
        std::string s1 = "silent";
        std::string s2 = argv[6];
        if (s1 == s2) {
            *silent = true;
        }
    }

    if (!*silent) {
        std::cout << "continue? (y/n):";
        char c;
        std::cin >> c;
        if (c != 'y') {
            wrong_N:
            std::cout << "Enter new N (must be even ans >= 6):";
            int N_usr;
            std::cin >> N_usr;
            if (N_usr >= 6 && N_usr % 2 == 0) {
                *N = N_usr;
            } else {
                goto wrong_N;
            }
            wrong_JRANGE:
            std::cout << "Enter new J_START (J1/J2):";
            double J_START_usr;
            std::cin >> J_START_usr;
            std::cout << "Enter new J_END (J1/J2):";
            double J_END_usr;
            std::cin >> J_END_usr;
            std::cout << "Enter new J_COUNT (number of data-points):";
            int J_COUNT_usr;
            std::cin >> J_COUNT_usr;
            if (J_START_usr <= J_END_usr && J_COUNT_usr >= 1) {
                *J_START = J_START_usr;
                *J_END = J_END_usr;
                *J_COUNT = J_COUNT_usr;
            } else {
                goto wrong_JRANGE;
            }
        }
    }
}
