#include "spinInversion.h"

/////////////////////////////// spin inversion ///////////////////////////////

namespace spinInversion {

    double get_gk(int k, int N) {
        if (k == 0 || k == N/4) {
            return 2.0;
        } else {
            return 1.0;
        }
    }

    int getClass_set_m_n(int &m, int &n, const int mp, const int mz, const int mpz) {

        if (mp == -1 && mz == -1 && mpz == -1) {m = -1; n = -1; return 1;}
        else if (mp != -1 && mz == -1 && mpz == -1) {m = mp; n = -1; return 2;}
        else if (mp == -1 && mz != -1 && mpz == -1) {m = mz; n = -1; return 3;}
        else if (mp == -1 && mz == -1 && mpz != -1) {m = mpz; n = -1; return 4;}
        else {m = mp; n = mz; return 5;}
        //else if (mp != -1 && mz != -1 && mpz == mp + mz) {m = mp; n = mz; return 5;}
        //else {m = -1; n = -1; return -1;}

    }

    double getNa(const int m, const int n, const unsigned int R, const int sigma, const int p, const int z, const int k, const int c, const int N) {

        double k_moment = 4.0 * PI * (double) k / (double) N;
        double Na = 2.0 * (double) pow(N, 2) / get_gk(k, N) / (double) R;

        switch (c) {
            case 1:
                Na *= 1;
                break;
            case 2:
                Na *= 1.0 + (double) sigma * (double) p * std::cos(k_moment * (double) m);
                break;
            case 3:
                Na *= 1.0 + (double) z * std::cos(k_moment * (double) m);
                break;
            case 4:
                Na *= 1.0 + (double) sigma * (double) p * (double) z * std::cos(k_moment * (double) m);
                break;
            case 5:
                Na *= ( 1.0 + (double) sigma * (double) p * std::cos(k_moment * (double) m) )
                        * ( 1.0 + (double) z * std::cos(k_moment * (double) n) );
                break;
            default:
                break;
        }

        //std::cout << Na << "\n";

        return Na;

    }

    double helement(const int a, const int b, const int l, const int q, const int g, const int k, const int p, const int z, const std::vector<int> &R_vals,
                    const std::vector<int> &m_vals, const std::vector<int> &n_vals, const std::vector<int> &c_vals, int N) {

        // state b
        int sigma_a = R_vals.at(a) / abs(R_vals.at(a));
        int m_a = m_vals.at(a);
        int n_a = n_vals.at(a);
        int c_a = c_vals.at(a);
        unsigned int R_a = std::abs(R_vals.at(a));
        double Na = getNa(m_a, n_a, R_a, sigma_a, p, z, k, c_a, N);
        //std::cout << "sigma_a: " << sigma_a << ", m_a: " << m_a << ", n_a: " << n_a << ", c_a: " << c_a << ", R_a: " << R_a << ", Na: " << Na << "\n";
        // state b
        int sigma_b = R_vals.at(b) / abs(R_vals.at(b));
        int m_b = m_vals.at(b);
        int n_b = n_vals.at(b);
        int c_b = c_vals.at(b);
        unsigned int R_b = std::abs(R_vals.at(b));
        double Nb = getNa(m_b, n_b, R_b, sigma_b, p, z, k, c_b, N);
        //std::cout << "sigma_b: " << sigma_b << ", m_b: " << m_b << ", n_b: " << n_b << ", c_b: " << c_b << ", R_b: " << R_b << ", Nb: " << Nb << "\n";

        double k_moment = 4.0 * PI * (double) k / (double) N;
        double val = 0.5 * (double) pow(sigma_a * p, q) * (double) std::pow(z, g) * sqrt(Nb / Na);
        //std::cout << "Na: " << Na << ", Nb: " << Na << ", val: " << val << "\n";
        //std::cout << "val: " << val << "\n";

        if (sigma_a == sigma_b) {
            switch (c_b) {
                case 1:
                case 3:
                    val *= std::cos(k_moment * (double) l);
                    break;
                case 2:
                case 5:
                    val *= ( std::cos(k_moment * (double) l) + (double) sigma_a * (double) p * std::cos(k_moment * (double) (l-m_b)) )
                            / ( 1.0 + (double) sigma_a * (double) p * std::cos(k_moment * (double) m_b) );
                    break;
                case 4:
                    val *= ( std::cos(k_moment * (double) l) + (double) sigma_a * (double) p * (double) z * std::cos(k_moment * (double) (l-m_b)) )
                          / ( 1.0 + (double) sigma_a * (double) p * (double) z * std::cos(k_moment * (double) m_b) );
                    break;
                default:
                    break;
            }
        } else {
            switch (c_b) {
                case 1:
                case 3:
                    val *= - (double) sigma_a * std::sin(k_moment * (double) l);
                    break;
                case 2:
                case 5:
                    val *= ( - (double) sigma_a * std::sin(k_moment * (double) l) + (double) p * std::sin(k_moment * (double) (l-m_b)) )
                           / ( 1.0 - (double) sigma_a * (double) p * std::cos(k_moment * (double) m_b) );
                    break;
                case 4:
                    val *= ( - (double) sigma_a * std::sin(k_moment * (double) l) + (double) p * (double) z * std::sin(k_moment * (double) (l-m_b)) )
                           / ( 1.0 - (double) sigma_a * (double) p * (double) z * std::cos(k_moment * (double) m_b) );
                    break;
                default:
                    break;
            }
        }
        if (std::abs(val) < epsilon) {
            val = 0.0;
        }

        //std::cout << "helement returned: " << val << "\n\n";

        return val;

    }

    void fillHamiltonSIBlock(const double &J1, const double &J2, int k, int p, int z, const std::vector<int> &states,
                                 const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                                 const std::vector<int> &c_vals, Eigen::MatrixXd &hamiltonBlock, const int &N) {

        int statesCount = (int) states.size();

        for (int a; a < statesCount; a++) {
            int state_n = 1;
            if (a > 0 && states.at(a - 1) == states.at(a)) {
                continue;
            } else if (a < states.size() - 1 && states.at(a) == states.at(a + 1)) {
                state_n = 2;
            }

            int s = states.at(a);

            for (int n = 0; n < N/2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J1;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -=  0.25 * J1;
                    }
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int r = 0, l = 0, q = 0, g = 0;
                    representative(d, r, l, q, g, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J1 * helement(i, j, l, q, g, k, p, z, R_vals, m_vals, n_vals, c_vals, N);
                            }
                        }
                    }
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J2;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -= 0.25 * J2;
                    }
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    int r = 0, l = 0, q = 0, g = 0;
                    representative(d, r, l, q, g, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J2 * helement(i, j, l, q, g, k, p, z, R_vals, m_vals, n_vals, c_vals, N);
                            }
                        }
                    }
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) += 0.25 * J2;
                    }
                } else {
                    for (int i = a; i < a + state_n; i++) {
                        hamiltonBlock(i, i) -= 0.25 * J2;
                    }
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    int r = 0, l = 0, q = 0, g = 0;
                    representative(d, r, l, q, g, N);
                    int b = findState(states, r);
                    if (b >= 0) {
                        int m = 1;
                        if (b > 0 && states.at(b) == states.at(b - 1)) {
                            m = 2, b -= 1;
                        } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                            m = 2;
                        }
                        for (int j = b; j < b + m; j++) {
                            for (int i = a; i < a + state_n; i++) {
                                hamiltonBlock(i, j) += J2 * helement(i, j, l, q, g, k, p, z, R_vals, m_vals, n_vals, c_vals, N);
                            }
                        }
                    }
                }
            }
        }

    }

    void SIBlockSolver(const double &J1, const double &J2, int k, int p, int z, const std::vector<int> &states,
                           const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                           const std::vector<int> &c_vals, std::vector<double> *eiVals, std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N) {

        const int statesCount = (int) states.size();
        Eigen::MatrixXd hamiltonBlock = Eigen::MatrixXd::Zero(statesCount, statesCount);
        fillHamiltonSIBlock(J1, J2, k, p, z, states, R_vals, m_vals, n_vals, c_vals, hamiltonBlock, N);

    #if defined(showMatrix) || defined(saveMatrix)
        matrixBlocks->push_back(hamiltonBlock);
    #endif

        Eigen::MatrixXd hamiltonBlockTransposed = hamiltonBlock.transpose();
        if (!hamiltonBlock.isApprox(hamiltonBlockTransposed)) {
            std::cout << "\nALARM: k: " << k << ", p: " << p << ", z: " << z << "\n";
        } else {
            //std::cout << "all good\n";
        }

        //std::cout << "calculating eigenvalues...\n";
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonBlock);
        const Eigen::VectorXd &H1EiVal = solver.eigenvalues();
        for (double ev: H1EiVal) {
            eiVals->push_back(ev);
        }

    }

    void getEiVals(const double &J1, const double &J2, std::vector<double> *eiVals,
                   std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N, const int &SIZE) {

        std::vector<std::vector<int>> states_m(N+1);
        for (int s = 0; s < SIZE; s++) {
            states_m.at(bitSum(s, N)).push_back(s);
        }

        std::vector<int> states;
        std::vector<int> R_vals;
        std::vector<int> m_vals;
        std::vector<int> n_vals;
        std::vector<int> c_vals;

        int numberOfStates = 0;
        const int k_upper = N/4;

        for (int mag = 0; mag <= N; mag++) {
            //std::cout << "mag: " << mag << "\n";
            for (int k = 0; k <= k_upper; k++) {

                if (mag == N/2) {
                    for (int z : {-1, 1}) {
                        for (int p: {-1, 1}) {
                            //if (k != 0 && k != k_upper && p == -1) {continue;}
                            for (int s: states_m.at(mag)) {
                                for (int sigma: {-1, 1}) {
                                    int R, n, m, mp, mz, mpz;
                                    checkStateSI(s, R, mp, mz, mpz, k, N);
                                    int c = getClass_set_m_n(m, n, mp, mz, mpz);
                                    //std::cout << "case: " << c << "\n";
                                    //std::cout << c;
                                    if ((k == 0 || k == k_upper) && sigma == -1) {continue;}
                                    if (c == 2 || c == 4 || c == 5) {
                                        double Na = getNa(m, n, R, sigma, p, z, k, c, N);
                                        double Na_inv = getNa(m, n, R, -sigma, p, z, k, c, N);
                                        if (std::abs(Na) < epsilon) {R = -1;}
                                        if (sigma == -1 && std::abs(Na_inv) > epsilon) {R = -1;}
                                    } else if (c == 3) {
                                        double val = 1.0 + (double) z * std::cos(4 * PI * (double) k * (double) m / (double) N);
                                        if (std::abs(val) < epsilon) {R = -1;}
                                    }
                                    if (R > 0) {
                                        //std::cout << getNa(m, n, R, sigma, p, z, k, c, N) << "\n";
                                        states.push_back(s);
                                        R_vals.push_back(sigma * R);
                                        m_vals.push_back(m);
                                        n_vals.push_back(n);
                                        c_vals.push_back(c);
                                        numberOfStates++;
                                    }
                                }
                            }
                            if (!states.empty()) {
                                SIBlockSolver(J1, J2, k, p, z, states, R_vals, m_vals, n_vals, c_vals, eiVals, matrixBlocks, N);
                            }
                            states.clear();
                            R_vals.clear();
                            m_vals.clear();
                            n_vals.clear();
                            c_vals.clear();
                        }
                    }
                }
                else {
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
                                    states.push_back(s);
                                    R_vals.push_back(sigma * R);
                                    m_vals.push_back(m);
                                    numberOfStates++;
                                }
                            }
                        }
                        if (!states.empty()) {
                            parityStates::parityBlockSolver(J1, J2, k, p, states, R_vals, m_vals, eiVals, matrixBlocks, N);
                        }
                        states.clear();
                        R_vals.clear();
                        m_vals.clear();
                    }
                }


            }
            //std::cout << "\n";
        }

        //std::cout << "number of states: " << numberOfStates << "\n";

        // sort eigenvalues
        std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });

    #if defined(showMatrix) || defined(saveMatrix)
        int offset_blocks = 0;
        Eigen::MatrixXd H_SI_Block = Eigen::MatrixXd::Zero(SIZE, SIZE);
        for (const Eigen::MatrixXd& M : *matrixBlocks) {
            H_SI_Block.block(offset_blocks, offset_blocks, M.rows(), M.cols()) = M;
            offset_blocks += (int) M.rows();
        }
    #endif
    #ifdef showMatrix
        std::cout << H_SI_Block << "\n";
    #endif
    #ifdef saveMatrix
        saveMatrixToFile(H_SI_Block, "HamiltonSpinInversion.txt", "spin inversion Ansatz für N = " + std::to_string(N) +
                                                               "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), N);
    #endif
    #ifdef showEigenvalues
        std::cout << "eigenvalues:\n";
        for (double ev : *eiVals) {
            std::cout << ev << "\n";
        }
    #endif
    #ifdef saveEigenvalues
        saveEiVals("EigenvaluesSpinInversion.txt", "spin inversion Ansatz für N = " + std::to_string(N) +
                                                       "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *eiVals, N);
    #endif
    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "spin inversion:..." << std::endl;

        auto *SIEiVals = new std::vector<double>;
        auto *matrixSIBlocks = new std::vector<Eigen::MatrixXd>;

        getEiVals(J1, J2, SIEiVals, matrixSIBlocks, N, SIZE);

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        delete SIEiVals;
        delete matrixSIBlocks;

    }

    void SIBlockSolver_withSave(const double &J1, const double &J2, int k, int p, int z, const std::vector<int> &states,
                                const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                                const std::vector<int> &c_vals, std::vector<double> &eiVals, std::vector<Eigen::MatrixXd> &matrixBlocks, const int &N) {

        const int statesCount = (int) states.size();
        Eigen::MatrixXd hamiltonBlock = Eigen::MatrixXd::Zero(statesCount, statesCount);
//        std::cout << "filling hamilton\n";
        fillHamiltonSIBlock(J1, J2, k, p, z, states, R_vals, m_vals, n_vals, c_vals, hamiltonBlock, N);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonBlock);
        const Eigen::VectorXd &H1EiVal = solver.eigenvalues();
        for (double ev: H1EiVal) {
            eiVals.push_back(ev);
        }

        eiVals.shrink_to_fit();

//        std::cout << "saving U Block\n";
        matrixBlocks.push_back(solver.eigenvectors());

//        // sort eigenvalues
//        std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
//            return std::real(c1) < std::real(c2);
//        });

    }

    void getEiValsZeroBlock(const double &J1, const double &J2, std::vector<double> & eiVals, std::vector<Eigen::MatrixXd> &UBlocks,
                            std::vector<int> &allStates, const int &N, const int &SIZE) {

        std::vector<int> states_m;
        fillStates(&states_m, N/2, N, SIZE);

        std::vector<int> states;
        std::vector<int> R_vals;
        std::vector<int> m_vals;
        std::vector<int> n_vals;
        std::vector<int> c_vals;

        int numberOfStates = 0;

        const int k_upper = N / 4;

        for (int k = 0; k <= k_upper; k++) {
            for (int z: {-1, 1}) {
                for (int p: {-1, 1}) {
                    for (int s: states_m) {
                        for (int sigma: {-1, 1}) {
                            int R, n, m, mp, mz, mpz;
                            checkStateSI(s, R, mp, mz, mpz, k, N);
                            int c = getClass_set_m_n(m, n, mp, mz, mpz);
                            if ((k == 0 || k == k_upper) && sigma == -1) { continue; }
                            if (c == 2 || c == 4 || c == 5) {
                                double Na = getNa(m, n, R, sigma, p, z, k, c, N);
                                double Na_inv = getNa(m, n, R, -sigma, p, z, k, c, N);
                                if (std::abs(Na) < epsilon) { R = -1; }
                                if (sigma == -1 && std::abs(Na_inv) > epsilon) { R = -1; }
                            } else if (c == 3) {
                                double val = 1.0 + (double) z * std::cos(4 * PI * (double) k * (double) m / (double) N);
                                if (std::abs(val) < epsilon) { R = -1; }
                            }
                            if (R > 0) {
                                allStates.push_back(s);
                                states.push_back(s);
                                R_vals.push_back(sigma * R);
                                m_vals.push_back(m);
                                n_vals.push_back(n);
                                c_vals.push_back(c);
                                numberOfStates++;
                            }
                        }
                    }
                    if (!states.empty()) {
//                        std::cout << "solving block\n";
                        SIBlockSolver_withSave(J1, J2, k, p, z, states, R_vals, m_vals, n_vals, c_vals, eiVals, UBlocks, N);
                    }
//                    std::cout << "clearing vectors\n";
                    states.clear();
                    R_vals.clear();
                    m_vals.clear();
                    n_vals.clear();
                    c_vals.clear();
                }
            }
        }

        std::cout << "number of states: " << numberOfStates << "\n";

    }

    //////////////////////// help ////////////////////////
    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "susceptibility (spin inversion): calculating..." << std::endl;

        std::vector<int> states;
        //fillStates(&states, N/2, N, SIZE);

        std::vector<double> eiVals;

        std::vector<Eigen::MatrixXd> UBlocks;

//        std::cout << "getting eiVals\n";
        getEiValsZeroBlock(J1, J2, eiVals, UBlocks, states, N, SIZE);

//        for (double ev : eiVals) {
//            std::cout << ev << "\n";
//        }
//
//        for (int a : states) {
//            printBits(a, N);
//        }
//
//        for (const Eigen::MatrixXd& M : UBlocks) {
//            std::cout << M << std::endl;
//        }

        ///// susceptibility /////

        const int statesCount = (int) states.size();

//        std::cout << "statesCount: " << statesCount << "\n";

        auto *susceptibility_magnetization = new std::vector<std::tuple<double, double>>;

//        std::cout << "getting S2\n";
        Eigen::MatrixXd S2 = spinMatrix(N, states);

//        std::cout << S2 << std::endl;

//        std::cout << "defining U\n";
        Eigen::MatrixXd U = Eigen::MatrixXd::Zero(statesCount, statesCount);

//        std::cout << U << std::endl;
//
//        std::cout << "getting U\n";

        int offset_blocks = 0;
        for (const Eigen::MatrixXd& M : UBlocks) {
            U.block(offset_blocks, offset_blocks, M.rows(), M.cols()) = M;
            offset_blocks += (int) M.rows();
        }

//        std::cout << U << std::endl;

//        std::cout << "getting U_inv_S2_U\n";
        Eigen::MatrixXd U_inv_S2_U = Eigen::MatrixXd::Zero(statesCount, statesCount);
        U_inv_S2_U = U.adjoint() * S2 * U;

//        std::cout << U_inv_S2_U << std::endl;

//        std::cout << "getting susceptibility" << std::endl;

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            //current_beta = 1 / current_beta;
            susceptibility_magnetization->push_back({current, getSusceptibilityDegeneracy(current, U_inv_S2_U, eiVals, N)});
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        ///// save /////

        std::string filenameSusceptibility_X = "spininversion_susceptibility_J_const.txt"; // spininversion_susceptibility_J_const.txt / data_susceptibility_J_const
        std::string headerSusceptibility_X = "N: " + std::to_string(N) + "\n"
                                             + "T START: " + std::to_string(START) + "\n"
                                             + "T END: " + std::to_string(END) + "\n"
                                             + "data-points: " + std::to_string(COUNT) + "\n"
                                             + "calculation time: " + formatTime(elapsed_seconds);

        std::string headerWithJSusceptibility_X = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSusceptibility_X;
        saveOutData(filenameSusceptibility_X, headerWithJSusceptibility_X, "J1/J2", "specific heat in J2", *susceptibility_magnetization, N);

        std::cout << "\n";

    }

    void startSpecificHeat(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                           const double &END, const int &COUNT) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "specific heat (spin inversion): calculating..." << std::endl;

        std::vector<double> eiVals;
        std::vector<Eigen::MatrixXd> matrixBlocks;
//        auto *eiVals = new std::vector<double>;
//        auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;

        getEiVals(J1, J2, &eiVals, &matrixBlocks, N, SIZE);


        ///// specific /////

//        auto *specificHeat_momentum = new std::vector<std::tuple<double, double>>;
        std::vector<std::tuple<double, double>> specificHeat_momentum;

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            specificHeat_momentum.emplace_back(current, getSpecificHeat(current, eiVals, N));
        }

        ///// save /////

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        std::string filenameSpecificHeat_C = "data_specific_heat_J_const.txt"; // spininversion_specific_heat
        std::string headerSpecificHeat_C = "N: " + std::to_string(N) + "\n"
                                           + "T START: " + std::to_string(START) + "\n"
                                           + "T END: " + std::to_string(END) + "\n"
                                           + "data-points: " + std::to_string(COUNT) + "\n"
                                           + "calculation time: " + formatTime(elapsed_seconds);

        std::string headerWithJSpecificHeat_C = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSpecificHeat_C;
        saveOutData(filenameSpecificHeat_C, headerWithJSpecificHeat_C, "J1/J2", "specific heat in J2", specificHeat_momentum, N);

        std::cout << "\n";

    }

}
