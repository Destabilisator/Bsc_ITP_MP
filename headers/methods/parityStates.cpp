#include "parityStates.h"

/////////////////////////////// parity states (unfinished) ///////////////////////////////

namespace parityStates {

    double get_gk(int k, int N) {
        if (k == 0 || k == N/4) {
            return 2.0;
        } else {
            return 1.0;
        }
    }

    double getNa(int sigma, int m, int Ra, int k, int p, int N) {
        //std::cout << sigma << " " << m << " " << Ra  << " " << k  << " " << p << "\n";
        double Na = (double) pow(N, 2) * get_gk(k, N) / (double) abs(Ra);
        //std::cout <<  "Na: " << Na << ", m = " << m << ", ";
        if (m != -1) {
            Na *= 1.0 + (double) sigma * (double) p * std::cos(4 * PI * (double) k * (double) m / (double) N); /// 4 * PI * (double) k * (double) m / (double) N
            //Na *= 1 + sigma * p * std::cos(4 * PI * (double) k * (double) m / (double) N);
            //std::cout << "\t" << Na << "\n";
        }
        //std::cout << "Na = " << Na << "\n";
        return Na;
    }

    double helement(int a, int b, int l, int q, int k, int p, const std::vector<int> &R_vals, const std::vector<int> &m_vals, int N) {
//        if (a == b) {
//            std::cout << "helement same state, ";
//        }
        int sigma_a = R_vals.at(a) / abs(R_vals.at(a));
        int m_a = m_vals.at(a);
        int R_a = R_vals.at(a);
        int sigma_b = R_vals.at(b) / abs(R_vals.at(b));
        int m_b = m_vals.at(b);
        int R_b = R_vals.at(b);
        double Na = getNa(sigma_a, m_a, R_a, k, p, N);
        double Nb = getNa(sigma_b, m_b, R_b, k, p, N);

        double k_moment = 4.0 * PI * (double) k / (double) N;

        //std::cout << "\t" << k << " " << l << " " << sigma_a << " " << sigma_b<< " " << p << " " << m_b << "\n";

        double val = 0.5 * (double) pow(sigma_a * p, q) * sqrt(Nb / Na);
        //std::cout << val << "\n";
        //std::cout << "\t" << Na << " " << Nb << " " << sigma_a << " " << p << " " << q << " " << pow(sigma_a * p, q) << "\n";

        if (sigma_a == sigma_b) { // same sigma
            if (m_b != -1) {
                //std::cout << "same sigma, m != -1, changed from " << val;
                val *= (std::cos(k_moment * (double) l) + (double) sigma_a * (double) p * std::cos(k_moment * (double) (l-m_b)))
                        / (1.0 + (double) sigma_a * (double) p * std::cos(k_moment * (double) m_b));
                //std::cout << " to " << val << "\n";
            } else {
                //std::cout << "same sigma changed from " << val;
                val *= std::cos(k_moment * (double) l);
                //std::cout << " to " << val << "\n";
            }
        } else {
            if (m_b != -1) {
                //std::cout << "different sigma, m != -1, changed from " << val;
                val *= (- (double) sigma_a * std::sin(k_moment * (double) l) + (double) p * std::sin(k_moment * (double) (l-m_b)))
                        / (1.0 - (double) sigma_a * (double) p * std::cos(k_moment * (double) m_b));
                //std::cout << " to " << val << "\n";
            } else {
                //std::cout << "different sigma changed from " << val;
                val *= - (double) sigma_a * std::sin(k_moment * (double) l);
                //std::cout << " to " << val << "\n";
            }
        }
        if (std::abs(val) < 1e-12) {
            val = 0.0;
        }
//        if ((a == 0 && b == 4) || (a == 4 && b == 0) || (a == 0 && b == 15) || (a == 15 && b == 0)) {
//            std::cout << "l = " << l << ", m_b = " << m_b << ", q = " << q << ", k = " << k << ", p = " << p << "\n";
//            std::cout << "helement returned: " << val << " at (" << a << "," << b << ")\n\n";
//        }
        return val;
    }

    void fillHamiltonParityBlock(const double &J1, const double &J2, const int &k, const int &p, const std::vector<int> &states,
                                 const std::vector<int> &R_vals, const std::vector<int> &m_vals,
                                 //Eigen::MatrixXcd &hamiltonBlock, const int &N) {
                                 Eigen::MatrixXd &hamiltonBlock, const int &N) {

        int statesCount = (int) states.size();

        for (int a; a < statesCount; a++) {
            //std::cout << "a = " << a << "\n";
            int state_n = 1;
            if (a > 0 && states.at(a - 1) == states.at(a)) {
                continue;
            } else if (a < states.size() - 1 && states.at(a) == states.at(a + 1)) {
                state_n = 2;
            }

            int s = states.at(a);

            for (int n = 0; n < N/2; n++) {
                //std::cout << "n = " << n << "\n";
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    //hamiltonBlock(a, a) += (double) state_n * std::complex<double>(0.25 * J1, 0.0);
                    //hamiltonBlock(a, a) += (double) state_n * 0.25 * J1;
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J1;
                    }
                } else {
                    //hamiltonBlock(a, a) -= (double) state_n * std::complex<double>(0.25 * J1, 0.0);
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -=  0.25 * J1;
                    }
                    //hamiltonBlock(a, a) -= (double) state_n * 0.25 * J1;
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int r = 0, l = 0, q = 0;
                    representative(d, &r, &l, &q, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        //std::cout << "j-loop 0-2\n";
                        for (int j = b; j < b + m; j++) {
                            //std::cout << "i-loop 0-2\n";
                            for (int i = a; i < a + state_n; i++) {
                                //if (i == j) { continue;}
                                //hamiltonBlock(i, j) += J1 * (std::complex<double>) helement(i, j, l, q, k, p, R_vals, m_vals, N);
                                hamiltonBlock(i, j) += J1 * helement(i, j, l, q, k, p, R_vals, m_vals, N);
                            }
                        }
                    }
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    //hamiltonBlock(a, a) += (double) state_n * std::complex<double>(0.25 * J2, 0.0);
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J2;
                    }

                    //hamiltonBlock(a, a) += (double) state_n * 0.25 * J2;
                } else {
                    //hamiltonBlock(a, a) -= (double) state_n * std::complex<double>(0.25 * J2, 0.0);
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -= 0.25 * J2;
                    }
                    //hamiltonBlock(a, a) -= (double) state_n * 0.25 * J2;
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    int r = 0, l = 0, q = 0;
                    representative(d, &r, &l, &q, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        //std::cout << "j-loop 0-1\n";
                        for (int j = b; j < b + m; j++) {
                            //std::cout << "i-loop 0-1\n";
                            for (int i = a; i < a + state_n; i++) {
                                //if (i == j) { continue;}
                                //hamiltonBlock(i, j) += J2 * (std::complex<double>) helement(i, j, l, q, k, p, R_vals, m_vals, N);
                                hamiltonBlock(i, j) += J2 * helement(i, j, l, q, k, p, R_vals, m_vals, N);
                            }
                        }
                    }
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    //hamiltonBlock(a, a) += (double) state_n * std::complex<double>(0.25 * J2, 0.0);
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J2;
                    }
                    //hamiltonBlock(a, a) += (double) state_n * 0.25 * J2;
                } else {
                    //hamiltonBlock(a, a) -= (double) state_n * std::complex<double>(0.25 * J2, 0.0);
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -= 0.25 * J2;
                    }
                    //hamiltonBlock(a, a) -= (double) state_n * 0.25 * J2;
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    int r = 0, l = 0, q = 0;
                    representative(d, &r, &l, &q, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        //std::cout << "j-loop 1-2\n";
                        for (int j = b; j < b + m; j++) {
                            //std::cout << "i-loop 1-2\n";
                            for (int i = a; i < a + state_n; i++) {
                                //if (i == j) { continue;}
                                //hamiltonBlock(i, j) += J2 * (std::complex<double>) helement(i, j, l, q, k, p, R_vals, m_vals, N);
                                hamiltonBlock(i, j) += J2 * helement(i, j, l, q, k, p, R_vals, m_vals, N);
                            }
                        }
                    }
                }
            }
        }

    }

    void parityBlockSolver(const double &J1, const double &J2, const int &k, const int &p, const std::vector<int> &states,
                           //const std::vector<int> &R_vals, const std::vector<int> &m_vals, std::vector<std::complex<double>> *eiVals,
                           //std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N) {
                           const std::vector<int> &R_vals, const std::vector<int> &m_vals, std::vector<double> *eiVals,
                           std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N) {

        const int statesCount = (int) states.size();
        //Eigen::MatrixXd hamiltonBlock(statesCount, statesCount);
        Eigen::MatrixXd hamiltonBlock = Eigen::MatrixXd::Zero(statesCount, statesCount);
        fillHamiltonParityBlock(J1, J2, k, p, states, R_vals, m_vals, hamiltonBlock, N);

    #if defined(showMatrix) || defined(saveMatrix)
        matrixBlocks->push_back(hamiltonBlock);
    #endif

        Eigen::MatrixXd hamiltonBlockTransposed = hamiltonBlock.transpose();

//        if (!hamiltonBlock.isApprox(hamiltonBlockTransposed, 1e-10)) {
//            std::cout << "ALARM!!!!\n";
//            matrixBlocks->push_back(hamiltonBlock);
//        }

        //std::cout << "calculating eigenvalues...\n";
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonBlock);
        const Eigen::VectorXd &H1EiVal = solver.eigenvalues();
        //for (std::complex<double> ev: H1EiVal) {
        for (double ev: H1EiVal) {
            eiVals->push_back(ev);
        }

    }

//    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *eiVals,
//                   std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE) {
    void getEiVals(const double &J1, const double &J2, std::vector<double> *eiVals,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE) {

        std::vector<std::vector<int>> states_m(N+1);
        for (int s = 0; s < SIZE; s++) {
            states_m.at(bitSum(s, N)).push_back(s);
        }

        std::vector<int> states;
        std::vector<int> R_vals;
        std::vector<int> m_vals;

        int numberOfStates = 0;
        const int k_upper = N/4;

        for (int mag = 0; mag <= N; mag++) {
            for (int k = 0; k <= k_upper; k++) {
                /////
                //if (mag != N/2 || k != 1) { continue;}////////////////////////////////////////////////////////////////////////
                /////
                for (int p : {-1, 1}) {
                    //if (k != 0 && k != k_upper && p == -1) {continue;}
                    for (int s : states_m.at(mag)) {
                        for (int sigma : {-1, 1}) {
                            int R, m;
                            checkState(s, &R, &m, k, N);
                            if ((k == 0 || k == k_upper) && sigma == -1) {continue;}
                            if (m != -1) {
                                double val = (double) sigma * (double) p * std::cos(4 * PI * (double) k * (double) m / (double) N);
                                if (std::abs(1.0 + val) < 1e-10) {R = -1;}
                                if (sigma == -1 && abs(1.0 - val) > 1e-10) {R = -1;}
                            }
                            if (R > 0) {
                                //std::cout << "m = " << mag << ", k = " << k << ", p = " << p << ", sigma = " << sigma << ": ";
                                //printBits(s, N);
                                states.push_back(s);
                                R_vals.push_back(sigma * R);
                                m_vals.push_back(m);
                                numberOfStates++;
                            }
                        }
                    }
                    if (!states.empty()) {
                        //std::cout << "new Block\n";
                        parityBlockSolver(J1, J2, k, p, states, R_vals, m_vals, eiVals, matrixBlocks, N);
                    }
                    states.clear();
                    R_vals.clear();
                    m_vals.clear();
                }
            }
        }

        //std::cout << "Number of states: " << numberOfStates << "\n";

        // sort eigenvalues
        std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });

#if defined(showMatrix) || defined(saveMatrix)
        int offset_blocks = 0;
//        Eigen::MatrixXd H_parity_Block = Eigen::MatrixXd::Zero(numberOfStates, numberOfStates);
        Eigen::MatrixXd H_parity_Block = Eigen::MatrixXd::Zero(SIZE, SIZE);
        for (const Eigen::MatrixXd& M : *matrixBlocks) {
            H_parity_Block.block(offset_blocks, offset_blocks, M.rows(), M.cols()) = M;
            offset_blocks += (int) M.rows();
        }
#endif
#ifdef showMatrix
        std::cout << H_parity_Block << "\n";
#endif
#ifdef saveMatrix
        saveMatrixToFile(H_parity_Block, "HamiltonParityStates.txt", "parity states Ansatz für N = " + std::to_string(N) +
                                                               "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), N);
#endif
#ifdef showEigenvalues
        std::cout << "eigenvalues:\n";
        for (double ev : *eiVals) {
            std::cout << ev << "\n";
        }
#endif
#ifdef saveEigenvalues
        //saveComplexEiVals("EigenvaluesParityStates.txt", "parity states Ansatz für N = " + std::to_string(N) +
        saveEiVals("EigenvaluesParityStates.txt", "parity states Ansatz für N = " + std::to_string(N) +
                                                       "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *eiVals, N);
#endif
    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {
        const clock_t begin_time_PARITY = clock();

        std::cout << "\nparity states:..." << std::endl;

        auto *parityEiVals = new std::vector<double>;
        auto *matrixParityBlocks = new std::vector<Eigen::MatrixXd>;

        getEiVals(J1, J2, parityEiVals, matrixParityBlocks, N, SIZE);

        auto time_PARITY = float(clock () - begin_time_PARITY) /  CLOCKS_PER_SEC;
        std::cout << "calculations done; this took: " << time_PARITY << " seconds\n";
        delete parityEiVals;
        delete matrixParityBlocks;
    }
}