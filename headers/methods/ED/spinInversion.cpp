#include "spinInversion.h"

/////////////////////////////// spin inversion ///////////////////////////////

namespace ED::spinInversion {

    double get_gk(int k, int N) {
        if (k == 0 || k == N/4) {
            return 2.0;
        } else {
            return 1.0;
        }
    }

    int getClass_set_m_n(int &m, int &n, const int &mp, const int &mz, const int &mpz) {

        if (mp == -1 && mz == -1 && mpz == -1) {m = -1; n = -1; return 1;}
        else if (mp != -1 && mz == -1 && mpz == -1) {m = mp; n = -1; return 2;}
        else if (mp == -1 && mz != -1 && mpz == -1) {m = mz; n = -1; return 3;}
        else if (mp == -1 && mz == -1 && mpz != -1) {m = mpz; n = -1; return 4;}
        else {m = mp; n = mz; return 5;}

    }

    double getNa(const int &m, const int &n, const unsigned int &R, const int &sigma, const int &p, const int &z, const int &k, const int &c, const int &N) {

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

        return Na;

    }

    double helement(const int &a, const int &b, const int &l, const int &q, const int &g, const int &k, const int &p,
                    const int &z, const std::vector<int> &R_vals, const std::vector<int> &m_vals,
                    const std::vector<int> &n_vals, const std::vector<int> &c_vals, int N) {

        // state b
        int sigma_a = R_vals.at(a) / abs(R_vals.at(a));
        int m_a = m_vals.at(a);
        int n_a = n_vals.at(a);
        int c_a = c_vals.at(a);
        unsigned int R_a = std::abs(R_vals.at(a));
        double Na = getNa(m_a, n_a, R_a, sigma_a, p, z, k, c_a, N);
        // state b
        int sigma_b = R_vals.at(b) / abs(R_vals.at(b));
        int m_b = m_vals.at(b);
        int n_b = n_vals.at(b);
        int c_b = c_vals.at(b);
        unsigned int R_b = std::abs(R_vals.at(b));
        double Nb = getNa(m_b, n_b, R_b, sigma_b, p, z, k, c_b, N);

        double k_moment = 4.0 * PI * (double) k / (double) N;
        double val = 0.5 * (double) pow(sigma_a * p, q) * (double) std::pow(z, g) * sqrt(Nb / Na);

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
        if (std::abs(val) < EPSILON) {
            val = 0.0;
        }

        return val;

    }

    void fillHamiltonSIBlock(const double &J1, const double &J2, const double &h, int k, int p, int z, const std::vector<int> &states,
                             const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                             const std::vector<int> &c_vals, Eigen::MatrixXd &hamiltonBlock, const int &N) {

        int statesCount = (int) states.size();

        for (int a; a < statesCount; a++) {

            int s = states.at(a);

            // external magnetic field
//            std::cout << hamiltonBlock(a,a) << std::endl;
//            int mag = bitSum(s, N) - N/2;
//            hamiltonBlock(a,a) += (double) mag * h;
//            std::cout << hamiltonBlock(a,a) << std::endl << std::endl;

            int state_n = 1;
            if (a > 0 && states.at(a - 1) == states.at(a)) {
                continue;
            } else if (a < states.size() - 1 && states.at(a) == states.at(a + 1)) {
                state_n = 2;
            }

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

    void SIBlockSolver(const double &J1, const double &J2, const double &h, int k, int p, int z, const std::vector<int> &states,
                           const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                           const std::vector<int> &c_vals, std::vector<double> *eiVals, std::vector<Eigen::MatrixXd> *matrixBlocks, const int &N) {

        const int statesCount = (int) states.size();
        Eigen::MatrixXd hamiltonBlock = Eigen::MatrixXd::Zero(statesCount, statesCount);
        fillHamiltonSIBlock(J1, J2, h, k, p, z, states, R_vals, m_vals, n_vals, c_vals, hamiltonBlock, N);

    #if defined(showMatrix) || defined(saveMatrix)
        matrixBlocks->push_back(hamiltonBlock);
    #endif

        //std::cout << "calculating eigenvalues...\n";
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonBlock);
        const Eigen::VectorXd &H1EiVal = solver.eigenvalues();
        for (double ev: H1EiVal) {
            eiVals->push_back(ev);
        }

    }

    void getEiVals(const double &J1, const double &J2, const double &h, std::vector<double> *eiVals,
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
            for (int k = 0; k <= k_upper; k++) {

                if (mag == N/2) {
                    for (int z : {-1, 1}) {
                        for (int p: {-1, 1}) {
                            for (int s: states_m.at(mag)) {
                                for (int sigma: {-1, 1}) {
                                    int R, n, m, mp, mz, mpz;
                                    checkStateSI(s, R, mp, mz, mpz, k, N);
                                    int c = getClass_set_m_n(m, n, mp, mz, mpz);
                                    if ((k == 0 || k == k_upper) && sigma == -1) {continue;}
                                    if (c == 2 || c == 4 || c == 5) {
                                        double Na = getNa(m, n, R, sigma, p, z, k, c, N);
                                        double Na_inv = getNa(m, n, R, -sigma, p, z, k, c, N);
                                        if (std::abs(Na) < EPSILON) { R = -1;}
                                        if (sigma == -1 && std::abs(Na_inv) > EPSILON) { R = -1;}
                                    } else if (c == 3) {
                                        double val = 1.0 + (double) z * std::cos(4 * PI * (double) k * (double) m / (double) N);
                                        if (std::abs(val) < EPSILON) { R = -1;}
                                    }
                                    if (R > 0) {
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
                                SIBlockSolver(J1, J2, h, k, p, z, states, R_vals, m_vals, n_vals, c_vals, eiVals, matrixBlocks, N);
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
                            parityStates::parityBlockSolver(J1, J2, h, k, p, states, R_vals, m_vals, eiVals, matrixBlocks, N);
                        }
                        states.clear();
                        R_vals.clear();
                        m_vals.clear();
                    }
                }


            }

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

    void start(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "spin inversion:..." << std::endl;

        auto *SIEiVals = new std::vector<double>;
        auto *matrixSIBlocks = new std::vector<Eigen::MatrixXd>;

        getEiVals(J1, J2, h, SIEiVals, matrixSIBlocks, N, SIZE);

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        delete SIEiVals;
        delete matrixSIBlocks;

    }

    //////////////////////// susceptibility ////////////////////////

    Eigen::MatrixXd spinMatrixSI(const int &N, int k, int p, int z, const std::vector<int> &states,
                                 const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                                 const std::vector<int> &c_vals) {

        int size = (int) states.size();
//        std::cout << "size: " << size << "\n";
        Eigen::MatrixXd S2 = 0.75 * (double) N * Eigen::MatrixXd::Identity(size, size);
        for (int a = 0; a < size; a++) {

            int state_n = 1;
            if (a > 0 && states.at(a - 1) == states.at(a)) {
                continue;
            } else if (a < states.size() - 1 && states.at(a) == states.at(a + 1)) {
                state_n = 2;
            }

            int s = states.at(a);

            for (int j = 0; j < N; j++) {
                for (int i = 0; i < j; i++) {
                    if (((s >> i) & 1) == ((s >> j) & 1)) {
                        for (int _ = a; _ < a + state_n; _++) {
//                            std::cout << "writing to " << a << " " << a << std::endl;
                            S2(a, a) += 0.5;
                        }
                    } else {
                        for (int _ = a; _ < a + state_n; _++) {
//                            std::cout << "writing to " << a << " " << a << std::endl;
                            S2(a, a) -= 0.5;
                        }
                        int d = s ^ (1 << i) ^ (1 << j);
                        int r = 0, l = 0, q = 0, g = 0;
                        representative(d, r, l, q, g, N);
                        int b = findState(states, r);
                        if (b >= 0) {
//                        S2(k, b) = 1.0;
                            int m = 1;
                            if (b > 0 && states.at(b) == states.at(b - 1)) {
                                m = 2, b -= 1;
                            } else if (b < states.size() - 1 && states.at(b) == states.at(b + 1)) {
                                m = 2;
                            }
                            for (int jj = b; jj < b + m; jj++) {
                                for (int ii = a; ii < a + state_n; ii++) {
//                                    std::cout << "writing to " << ii << " " << jj << std::endl;
                                    S2(ii, jj) += 2.0 * helement(ii, jj, l, q, g, k, p, z, R_vals, m_vals, n_vals, c_vals, N);
                                }
                            }
                        }
                    }
                }
            }
        }
//        for (int i = 0; i < size; i++) {
//            std::cout << S2(i,i) << std::endl;
//        }
        //std::cout << S2 << std::endl;
        return S2;
    }

    void SIBlockSolver_withMatrix(const double &J1, const double &J2, const double &h, int k, int p, int z, const std::vector<int> &states,
                                  const std::vector<int> &R_vals, const std::vector<int> &m_vals, const std::vector<int> &n_vals,
                                  const std::vector<int> &c_vals, std::vector<std::vector<double>> &eiVals, std::vector<Eigen::MatrixXd> &matrixUBlocks,
                                  std::vector<Eigen::MatrixXd> &matrixS2Blocks, const int &N) {

        const int statesCount = (int) states.size();
        Eigen::MatrixXd hamiltonBlock = Eigen::MatrixXd::Zero(statesCount, statesCount);
        fillHamiltonSIBlock(J1, J2, h, k, p, z, states, R_vals, m_vals, n_vals, c_vals, hamiltonBlock, N);

//        std::cout << "hamiltonblock\n";
//        std::cout << hamiltonBlock << std::endl;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonBlock);
        const Eigen::VectorXd &H1EiVal = solver.eigenvalues();
//        std::cout << "egenvalues: \n";
        std::vector<double> eV;
        for (double ev: H1EiVal) {
            eV.push_back(ev);
//            std::cout << ev << "\n";
        }
        eV.shrink_to_fit();
        eiVals.push_back(eV);
        eiVals.shrink_to_fit();

        matrixUBlocks.push_back(solver.eigenvectors());

//        std::cout << "defining S2" << std::endl;
        Eigen::MatrixXd S2 = Eigen::MatrixXd::Zero(statesCount, statesCount);
//        std::cout << "filling S2" << std::endl;
        S2 = spinMatrixSI(N, k, p, z, states, R_vals, m_vals, n_vals, c_vals);
//        std::cout << "saving S2" << std::endl;
        matrixS2Blocks.push_back(S2);

//        std::cout << "S2:" << std::endl;
//        std::cout << S2 << std::endl;

//        // sort eigenvalues
//        std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
//            return std::real(c1) < std::real(c2);
//        });

    }

    void getEiValsZeroBlock(const double &J1, const double &J2, const double &h, std::vector<std::vector<double>> &eiVals, std::vector<Eigen::MatrixXd> &UBlocks,
                            std::vector<Eigen::MatrixXd> &S2Blocks, const int &N, const int &SIZE) {

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
                                if (std::abs(Na) < EPSILON) { R = -1; }
                                if (sigma == -1 && std::abs(Na_inv) > EPSILON) { R = -1; }
                            } else if (c == 3) {
                                double val = 1.0 + (double) z * std::cos(4 * PI * (double) k * (double) m / (double) N);
                                if (std::abs(val) < EPSILON) { R = -1; }
                            }
                            if (R > 0) {
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
//                        std::cout << "block\n";
                        SIBlockSolver_withMatrix(J1, J2, h, k, p, z, states, R_vals, m_vals, n_vals, c_vals, eiVals,
                                                 UBlocks, S2Blocks, N);
                    }
                    states.clear();
                    R_vals.clear();
                    m_vals.clear();
                    n_vals.clear();
                    c_vals.clear();
                }
            }
        }

        //std::cout << "number of states: " << numberOfStates << "\n";

    }

    void getEiValsMagBlock(const double &J1, const double &J2, const double &h, std::vector<double> *eiVals, const int &N, const int &SIZE, const int &mag) {

        std::vector<Eigen::MatrixXd> matrixBlocks;

        std::vector<int> states_m;
        fillStates(&states_m, mag, N, SIZE);

        std::vector<int> states;
        std::vector<int> R_vals;
        std::vector<int> m_vals;
        std::vector<int> n_vals;
        std::vector<int> c_vals;

        int numberOfStates = 0;
        const int k_upper = N/4;

        for (int k = 0; k <= k_upper; k++) {

            if (mag == N/2) {
                for (int z : {-1, 1}) {
                    for (int p: {-1, 1}) {
                        for (int s: states_m) {
                            for (int sigma: {-1, 1}) {
                                int R, n, m, mp, mz, mpz;
                                checkStateSI(s, R, mp, mz, mpz, k, N);
                                int c = getClass_set_m_n(m, n, mp, mz, mpz);
                                if ((k == 0 || k == k_upper) && sigma == -1) {continue;}
                                if (c == 2 || c == 4 || c == 5) {
                                    double Na = getNa(m, n, R, sigma, p, z, k, c, N);
                                    double Na_inv = getNa(m, n, R, -sigma, p, z, k, c, N);
                                    if (std::abs(Na) < EPSILON) { R = -1;}
                                    if (sigma == -1 && std::abs(Na_inv) > EPSILON) { R = -1;}
                                } else if (c == 3) {
                                    double val = 1.0 + (double) z * std::cos(4 * PI * (double) k * (double) m / (double) N);
                                    if (std::abs(val) < EPSILON) { R = -1;}
                                }
                                if (R > 0) {
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
                            SIBlockSolver(J1, J2, h, k, p, z, states, R_vals, m_vals, n_vals, c_vals, eiVals, &matrixBlocks, N);
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
                    for (int s : states_m) {
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
                        parityStates::parityBlockSolver(J1, J2, k, p, states, R_vals, m_vals, eiVals, &matrixBlocks, N);
                    }
                    states.clear();
                    R_vals.clear();
                    m_vals.clear();
                }
            }


        }

    }

    ///// help /////
    void startSusceptibility(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
                             const double &END, const int &COUNT) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "susceptibility (spin inversion): calculating..." << std::endl;

        std::vector<std::vector<double>> eiVals;

        std::vector<Eigen::MatrixXd> UBlocks;
        std::vector<Eigen::MatrixXd> S2Blocks;

        std::cout << "getting eiVals\n";
        getEiValsZeroBlock(J1, J2, 0.0, eiVals, UBlocks, S2Blocks, N, SIZE);

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

        //auto *susceptibility_magnetization = new std::vector<std::tuple<double, double>>;
        std::vector<std::tuple<double, double>> susceptibility_magnetization;

        std::vector<Eigen::MatrixXd> U_inv_S2_U;

//        std::cout << std::endl << std::endl;

        int matrixCount = (int) UBlocks.size();
        for (int i = 0; i < matrixCount; i++) {
            int sz = (int) UBlocks.at(i).rows();
            Eigen::MatrixXd M = Eigen::MatrixXd::Zero(sz, sz);
            M = UBlocks.at(i).transpose() * S2Blocks.at(i) * UBlocks.at(i);
            U_inv_S2_U.push_back(M);
//            std::cout << M << std::endl;
//            std::cout << M.rows() << std::endl;
        }

//        std::cout << std::endl << std::endl;
//        for (std::vector<double> eV : eiVals) {
//            std::cout << eV.size() << "\n";
//            std::cout << "new ev block\n";
//            for (double ev : eV) {
//                std::cout << ev << "\n";
//            }
//        }

        for (int i = 0; i <= 0; i++) { // for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            //current_beta = 1 / current;
            susceptibility_magnetization.emplace_back(current, getSusceptibilityDegeneracy(current, U_inv_S2_U, eiVals, N));
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

//        for (std::tuple<double, double> t : susceptibility_magnetization) {
//            std::cout << std::get<0>(t) << " " << std::get<1>(t) << "\n";
//        }

        ///// save /////

        std::string filenameSusceptibility_X = "spininversion_susceptibility_J_const.txt"; // spininversion_susceptibility_J_const.txt / data_susceptibility_J_const
        std::string headerSusceptibility_X = "N: " + std::to_string(N) + "\n"
                                             + "T START: " + std::to_string(START) + "\n"
                                             + "T END: " + std::to_string(END) + "\n"
                                             + "data-points: " + std::to_string(COUNT) + "\n"
                                             + "calculation time: " + formatTime(elapsed_seconds);

        std::string headerWithJSusceptibility_X = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSusceptibility_X;
        saveOutData(filenameSusceptibility_X, headerWithJSusceptibility_X, "J1/J2", "specific heat in J2", susceptibility_magnetization, N);
//        std::cout << "\n";

    }

    //////////////////////// specific heat ////////////////////////
    void startSpecificHeat(const double &J1, const double &J2, const double &h, const int &N, const int &SIZE, const double &START,
                           const double &END, const int &COUNT) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "specific heat (spin inversion): calculating..." << std::endl;

        std::vector<double> eiVals;
        std::vector<Eigen::MatrixXd> matrixBlocks;

        getEiVals(J1, J2, h, &eiVals, &matrixBlocks, N, SIZE);


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

        std::string filenameSpecificHeat_C = "data_specific_heat_J_const.txt"; // data_specific_heat_J_const_SI / data_specific_heat_J_const
        std::string headerSpecificHeat_C = "N: " + std::to_string(N) + "\n"
                                           + "h: " + std::to_string(h) + "\n"
                                           + "T START: " + std::to_string(START) + "\n"
                                           + "T END: " + std::to_string(END) + "\n"
                                           + "data-points: " + std::to_string(COUNT) + "\n"
                                           + "calculation time: " + formatTime(elapsed_seconds);

        std::string headerWithJSpecificHeat_C = "J1/J2: " + std::to_string(J1/J2) +"\n" + headerSpecificHeat_C;
        saveOutData(filenameSpecificHeat_C, headerWithJSpecificHeat_C, "J1/J2", "specific heat in J2", specificHeat_momentum, N);
//        std::cout << "\n";

    }

}
