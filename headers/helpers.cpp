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

//int invertBits(int s, int N) {
//    int t = 0;
//    for (int n = 0; n < N; n++) {
//        t |= ( (s >> n) & 1 ^ 1) << n;
//    }
//    return t;
//}

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
    for (int i = 0; i <= *r; i++) { // for (int i = 0; i < *r; i++) {
        if (t < s) {*r = -1; return;}
        else if (t == s) {*m = i; return;}
        //t = translateLeft(t, 1, N);
        t = translateLeft(t, 2, N);
    }
}

void checkStateSI(const int &s, int &r, int &mp, int &mz, int &mpz, const int &k, const int &N) {
    int t = s; r = -1;
    // no Z or P
    for (int i = 1; i <= N/2; i++) {
        t = translateLeft(t, 2, N);
        if (t < s) {return;}
        else if (t == s) {
            if (k % (int) ((double) N / (double) i / 2.0) != 0) {return;}
            else {r = i; break;}
        }
    }
    // only P
    t = reflectBits(translateLeft(s, 1, N), N); mp = -1;
    for (int i = 0; i <= r; i++) {
        if (t < s) {r = -1; return;}
        else if (t == s) {mp = i; break;}
        t = translateLeft(t, 2, N);
    }
    // only Z
    t = invertBits(s, N); mz = -1;
    for (int i = 0; i <= r; i++) {
        if (t < s) {r = -1; return;}
        else if (t == s) {mz = i; break;}
        t = translateLeft(t, 2, N);
    }
    // Z and P
    t = reflectBits(translateLeft(invertBits(s, N), 1, N), N); mpz = -1;
    for (int i = 0; i <= r; i++) {
        if (t < s) {r = -1; return;}
        else if (t == s) {mpz = i; break;}
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
    for (int i = 1; i <= N/2; i++) {
        t = translateLeft(t, 2, N);
        if (t < *r) {*r = t; *l = i;}
    }
    t = reflectBits(translateLeft(s, 1, N), N);
//    t = reflectBits(translateLeft(t, 1, N), N);
    *q = 0;
    for (int i = 0; i <= N/2; i++) {
        if (t < *r) {*r = t; *l = i; *q = 1;}
        t = translateLeft(t, 2, N);
    }
}

void representative(const int &s, int &r, int &l, int &q, int &g, const int &N) {
    int t = s; r = s; l = 0; q = 0; g = 0;
    // no Z or P
    for (int i = 1; i <= N/2; i++) {
        t = translateLeft(t, 2, N);
        if (t < r) {r = t; l = i; q = 0; g = 0;}
    }
    // only P
    t = reflectBits(translateLeft(s, 1, N), N);
    for (int i = 0; i <= N/2; i++) {
        if (t < r) {r = t; l = i; q = 1; g = 0;}
        t = translateLeft(t, 2, N);
    }
    // only Z
    t = invertBits(s, N);
    for (int i = 0; i <= N/2; i++) {
        if (t < r) {r = t; l = i; q = 0; g = 1;}
        t = translateLeft(t, 2, N);
    }
    // Z and P
    t = reflectBits(translateLeft(invertBits(s, N), 1, N), N);
    for (int i = 0; i <= N/2; i++) {
        if (t < r) {r = t; l = i; q = 1; g = 1;}
        t = translateLeft(t, 2, N);
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
                file << hamilton[i][j] << "\t" ;
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

void saveOutData(const std::string &filename, const std::string &header, const std::string &x_label,
                 const std::string &y_label, const std::vector<std::tuple<double, double, int, int>> &outData, const int &N) {

    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
        file << x_label << "\t" << y_label << "\n";
        for (std::tuple<double, double, int, int> data : outData) {
            file << std::get<0>(data) << "\t" << std::get<1>(data)
                    << "\t" << std::get<2>(data)  << "\t" << std::get<3>(data) << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void saveOutData(const std::string &filename, const std::string &header, const std::string &x_label,
                 const std::string &y_label, const std::vector<std::tuple<double, double, int, int, int, int>> &outData, const int &N) {

    std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << "\n";
    std::ofstream file;
    try {
        file.open("./results/" + std::to_string(N) + "_" + filename);
        file << header << "\n\n";
        file << x_label << "\t" << y_label << "\n";
        for (std::tuple<double, double, int, int, int, int> data : outData) {
            #ifdef fixedPrecision
                file << std::fixed << std::setprecision(5) << std::get<0>(data) << "\t" << std::get<1>(data);
            #else
                file << std::get<0>(data) << "\t" << std::get<1>(data);
            #endif
            file << "\t" << std::get<2>(data)  << "\t" << std::get<3>(data)
                 << "\t" << std::get<4>(data)  << "\t" << std::get<5>(data)<< "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void save3DPlotDataC(const double &J, const int &N, const std::vector<std::tuple<double, double>>& C_func_T) {
    std::ofstream file;
    try {
        file.open("./results/3DData/" + std::to_string(N) + "/C/" + std::to_string(J) + ".txt");
        for (std::tuple<double, double> data : C_func_T) {
            file << std::get<0>(data) << "\t" << std::get<1>(data) << "\n";
        }
    } catch (...) {
        file.close();
        std::cout << "failed to save to file\n";
    }
    file.close();
}

void save3DPlotDataX(const double &J, const int &N, const std::vector<std::tuple<double, double>>& X_func_T) {
    std::ofstream file;
    try {
        file.open("./results/3DData/" + std::to_string(N) + "/X/" + std::to_string(J) + ".txt");
        for (std::tuple<double, double> data : X_func_T) {
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

Eigen::MatrixXcd spinMatrixMomentum(const int &N, const int &k, const std::vector<int> &states, const std::vector<int> &R_vals) {

    const int size = (int) states.size();
    Eigen::MatrixXcd S2 = 0.75 * (double) N * Eigen::MatrixXd::Identity(size, size);
    for (int a = 0; a < size; a++) {
        int s = states.at(a);
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < j; i++) {
                if (((s >> i) & 1) == ((s >> j) & 1)) {
                    S2(a, a) += 0.5;
                } else {
                    S2(a, a) -= 0.5;
                    int d = s ^ (1 << i) ^ (1 << j);
                    int r = 0, l = 0;
                    representative(d, &r, &l, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                        S2(a,b) += (std::complex<double>) 1.0 * sqrt((double) R_vals.at(a) / (double) R_vals.at(b)) * std::exp(numC);
                    }
                }
            }
        }
    }
    //std::cout << S2 << std::endl;
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

double getSpecificHeat(const double &temp, const std::vector<double>& eiVals, const int &N) {
    double Z_sum = 0.0, expectation_H = 0.0, expectation_H_2 = 0.0;
    for (double ev : eiVals) {
        Z_sum += std::exp(-1.0 / temp * ev);
        expectation_H += std::exp(-1.0 / temp * ev) * ev;
        expectation_H_2 += std::exp(-1.0 / temp * ev) * ev * ev;
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

double getSusceptibilityDegeneracy(const double &temp, const Eigen::MatrixXd &M, const std::vector<double>& eiVals, const int &N) {
    double Z_sum = 0.0, expectation_mz_2 = 0.0;
    for (int i = 0; i < eiVals.size(); i++) {
        double ev_real = eiVals.at(i);
        double S_elem = M(i, i);
        double S = - 0.5 + std::sqrt(0.25 + S_elem);
        Z_sum += std::exp(-1.0 / temp * ev_real) * (2.0 * S + 1);
        expectation_mz_2 += std::exp(-1.0 / temp * ev_real) * S_elem * (2.0 * S + 1);
    }
    expectation_mz_2 /= Z_sum;
    expectation_mz_2 /= 3.0;
    return 1.0 / temp * expectation_mz_2 / N;
}

double getSusceptibilityDegeneracy(const double &temp, const std::vector<Eigen::MatrixXd> &M_list, const std::vector<std::vector<double>> &eiVal_list, const int &N) {
    double Z_sum = 0.0, expectation_mz_2 = 0.0;
    for (int cur = 0; cur < M_list.size(); cur ++) {
//        std::cout << "new block:\n";
        const Eigen::MatrixXd &M = M_list.at(cur);
        std::vector<double> eiVals = eiVal_list.at(cur);
        for (int i = 0; i < eiVals.size(); i++) {
            double ev_real = eiVals.at(i);
            double S_elem = M(i, i);//std::max(-0.25, M(i, i));
            double S = - 0.5 + std::sqrt(0.25 + S_elem);
//            std::cout << "S_elem: " << S_elem;
//            std::cout << ", S: " << S << std::endl;
            Z_sum += std::exp(-1.0 / temp * ev_real) * (2.0 * S + 1);
            expectation_mz_2 += std::exp(-1.0 / temp * ev_real) * S_elem * (2.0 * S + 1);
            //std::cout << ev_real << " " << S_elem << " " << S << ": " << Z_sum << " " << expectation_mz_2 << "\n";
        }
    }
    expectation_mz_2 /= Z_sum;
    expectation_mz_2 /= 3.0;
    //std::cout << 1.0 / temp * expectation_mz_2 / N << "\n\n\n";
    return 1.0 / temp * expectation_mz_2 / N;
}

double getSusceptibilityDegeneracy(const double &temp, const std::vector<Eigen::MatrixXcd> &M_list, const std::vector<std::vector<std::complex<double>>> &eiVal_list, const int &N) {
    double Z_sum = 0.0, expectation_mz_2 = 0.0;
    for (int cur = 0; cur < M_list.size(); cur ++) {
        const Eigen::MatrixXcd &M = M_list.at(cur);
        std::vector<std::complex<double>> eiVals = eiVal_list.at(cur);
        for (int i = 0; i < eiVals.size(); i++) {
            double ev_real = std::real(eiVals.at(i));
            double S_elem = std::real(M(i, i));
            double S = - 0.5 + std::sqrt(0.25 + S_elem);
            if (S < epsilon) {S = 0.0;}
//            std::cout << "S_elem: " << S_elem;
//            std::cout << ", S: " << S << std::endl;
            Z_sum += std::exp(-1.0 / temp * ev_real) * (2.0 * S + 1);
            expectation_mz_2 += std::exp(-1.0 / temp * ev_real) * S_elem * (2.0 * S + 1);
        }
    }
    expectation_mz_2 /= Z_sum;
    expectation_mz_2 /= 3.0;
    return 1.0 / temp * expectation_mz_2 / N;
}

/////////////////////////////// others ///////////////////////////////

// [executable] N J_START J_END J_COUNT CORES SILENT
void validateInput(int &argc, char* argv[], const unsigned int &cpu_cnt, int &N, int &SIZE, double &J_START, double &J_END,
                   int &J_COUNT, double &T_START, double &T_END, int &T_COUNT, bool &silent, int &cores, bool &plotsIn3D,
                   bool skipSilent, const double &J1, const double &J2, bool &noX) {

    if (argc >= 2) {
        std::string DDD = "3D";
        if (argv[1] != DDD) {
            goto no3D;
        }
        plotsIn3D = true;
        if ( (std::stoi(argv[2]) % 2 == 0) && (std::stoi(argv[2]) >= 6) && (std::stoi(argv[2]) <= 32) ) {
            N = std::stoi(argv[2]);
        } else {
            std::cout << "invalid chain size, must be even and at least 6, defaulting to " << N << "\n";
        }
        SIZE = (int) pow(2, N);
        std::cout << "N: " << N << "; size: " << SIZE << std::endl;
        if (std::stod(argv[3]) > std::stod(argv[4]) || std::stoi(argv[5]) < 1) {
            std::cout << "range invalid, defaulting...\n";
        } else {
            J_START = std::stod(argv[3]); J_END = std::stod(argv[4]); J_COUNT = std::stoi(argv[5]);
        }
        std::cout << "J_START = " << J_START << ", J_END = " << J_END << " and J_COUNT = " << J_COUNT << " from args\n";

        if (std::stod(argv[6]) > std::stod(argv[7]) || std::stoi(argv[8]) < 1) {
            std::cout << "range invalid, defaulting...\n";
        } else {
            T_START = std::stod(argv[6]); T_END = std::stod(argv[7]); T_COUNT = std::stoi(argv[8]);
        }
        std::cout << "T_START = " << T_START << ", T_END = " << T_END << " and T_COUNT = " << T_COUNT << " from args\n";

        if (argc >= 9) {
            int crs = std::stoi(argv[9]);
            if (crs > 0 && crs <= cpu_cnt) {
                cores = crs;
                std::cout << "using " << cores << " cores\n";
            } else {
                std::cout << "defaulting to using all (" << cores << ") cores\n";
            }
        }

        if (argc >= 10) {
            std::string s1 = "noX";
            std::string s2 = argv[10];
            std::cout << s2 << "\n";
            if (s1 == s2) {
                noX = true;
            }
        }

        if (skipSilent) {
            return;
        }

        if (argc >= 11) {
            std::string s1 = "silent";
            std::string s2 = argv[11];
            std::cout << s2 << "\n";
            if (s1 == s2) {
                silent = true;
            }
        }

        goto silentEntry;
    }

    no3D:

    if (argc >= 2) {
        if ( (std::stoi(argv[1]) % 2 == 0) && (std::stoi(argv[1]) >= 6) && (std::stoi(argv[1]) <= 32) ) {
            N = std::stoi(argv[1]);
        } else {
            std::cout << "invalid chain size, must be even and at least 6, defaulting to " << N << "\n";
        }
    }
    SIZE = (int) pow(2, N);
    std::cout << "N: " << N << "; size: " << SIZE << std::endl;
    if (argc >= 5) {
        std::cout << "range given: ";
        if (std::stod(argv[2]) > std::stod(argv[3]) || std::stoi(argv[4]) < 1) {
            std::cout << "range invalid, defaulting...\n";
        } else {
            J_START = std::stod(argv[2]); J_END = std::stod(argv[3]); J_COUNT = std::stoi(argv[4]);
        }
        std::cout << "START = " << J_START << ", END = " << J_END << " and COUNT = " << J_COUNT << " from args\n";
    } else {
        std::cout << "no range given: ";
        std::cout << "START = " << J_START << ", END = " << J_END << " and COUNT = " << J_COUNT << " from default\n";
    }

    std::cout << "default J1 and J2 for beta plot: J1 = " << J1 << " and J2 = " << J2 << " (currently unchangeable)\n";

    if (argc >= 6) {
        int crs = std::stoi(argv[5]);
        if (crs > 0 && crs <= cpu_cnt) {
            cores = crs;
            std::cout << "using " << cores << " cores\n";
        } else {
            std::cout << "defaulting to using all (" << cores << ") cores\n";
        }
    }

    T_START = J_START;
    T_END = J_END;
    T_COUNT = J_COUNT;

    if (argc >= 7) {
        std::string s1 = "noX";
        std::string s2 = argv[6];
        if (s1 == s2) {
            noX = true;
        }
    }

    if (skipSilent) {
        return;
    }

    if (argc >= 8) {
        std::string s1 = "silent";
        std::string s2 = argv[7];
        if (s1 == s2) {
            silent = true;
        }
    }

    silentEntry:

    if (!silent) {
        std::cout << "continue? (y/n):";
        char c;
        std::cin >> c;
        if (c != 'y') {
            exit(1);
//            wrong_N:
//            std::cout << "Enter new N (must be even ans >= 6):";
//            int N_usr;
//            std::cin >> N_usr;
//            if (N_usr >= 6 && N_usr % 2 == 0) {
//                N = N_usr;
//            } else {
//                goto wrong_N;
//            }
//            wrong_JRANGE:
//            std::cout << "Enter new J_START (J1/J2):";
//            double J_START_usr;
//            std::cin >> J_START_usr;
//            std::cout << "Enter new J_END (J1/J2):";
//            double J_END_usr;
//            std::cin >> J_END_usr;
//            std::cout << "Enter new J_COUNT (number of data-points):";
//            int J_COUNT_usr;
//            std::cin >> J_COUNT_usr;
//            if (J_START_usr <= J_END_usr && J_COUNT_usr >= 1) {
//                J_START = J_START_usr;
//                J_END = J_END_usr;
//                J_COUNT = J_COUNT_usr;
//            } else {
//                goto wrong_JRANGE;
//            }
        }
    }
}

std::string formatTime(std::chrono::duration<double> elapsed_seconds) {
    double raw = elapsed_seconds.count();
    int hours = (int) raw / 3600;
    int minutes = (int) (raw - hours * 3600) / 60;
    double seconds = raw - hours * 3600 - minutes *  60;
    std::string ret;
    if (hours > 0) {ret += std::to_string(hours) + " hours " + std::to_string(minutes) + " minutes ";}
    else if (minutes > 0) {ret += std::to_string(minutes) + " minutes ";}
    ret += std::to_string(seconds) + " seconds";
    return ret;
}
