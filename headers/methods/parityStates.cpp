#include "parityStates.h"

/////////////////////////////// parity states (unfinished) ///////////////////////////////

namespace parityStates {
    void getEiVals(double J1, double J2, std::vector<std::complex<double>> *eiVals,
                   std::vector<Eigen::MatrixXcd> *matrixBlocks, const int &N, const int &SIZE) {

        int k_lower = -(N+2)/4+1;
        int k_upper = N/4;

        std::vector<std::vector<std::vector<int>>> vec(N+1, std::vector<std::vector<int>>(N));
        auto *states_mk = &vec;

        for (int s = 0; s < SIZE; s++) {
            int m = bitSum(s, N);
            for (int k = k_lower; k <= k_upper; k++) {
                int R = checkState(s, k, N);
                if (R >= 0) {
                    states_mk->at(m).at(k-k_lower).push_back(s);
                }
            }
        }

        auto *states = new std::vector<int>;
        auto *R_vals = new std::vector<int>;
        auto *m_vals = new std::vector<int>;

        for (int m = 0; m <= N; m++) {
            for (int k = k_lower; k <= k_upper; k++) {
                for (int s : states_mk->at(m).at(k-k_lower)) {
                    int R, m_cs;
                    checkState(s, &R, &m_cs, k, N);
                    for (int p : {-1, 1}) {
                        // p = +1 for k != 0, pi
                        if (k != 0 && k != k_upper && p == -1 ) {
                            continue;
                        }
                        for (int sigma : {-1, 1}) {
                            // sigma = +1 if k == 0 or k == k_upper (N/4)
                            if ((k == 0 || k == k_upper) && sigma == -1) {
                                continue;
                            }
                            if (m_cs != -1) {
                                std::complex<double> val  = (double) sigma * (double) p * std::cos(std::complex<double>(0, 4 * PI * (double) k * (double) m_cs / (double) N));
                                if (abs(std::complex<double>(1,0) + val) < 0.0001) {R = -1;}
                                if (sigma == -1 && abs(std::complex<double>(1,0) - val) > 0.0001) {R = -1;}
                            } if (R > 0) {
                                states->push_back(s);
                                R_vals->push_back(sigma * R);
                                m_vals->push_back(m_cs);
                            }

                        }
                        // ??????????
                    }
                }


                //parityBlockSolver(J1, J2, k, states->at(m).at(k-k_lower), R_vals->at(m).at(k-k_lower), eiVals, matrixBlocks);
            }
        }

        // sort eigenvalues
        std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });

#if defined(showMatrix) || defined(saveMatrix)
        int offset_blocks = 0;
    Eigen::MatrixXcd H_parity_Block = Eigen::MatrixXcd::Zero(SIZE, SIZE);
    for (const Eigen::MatrixXcd& M : *matrixBlocks) {
        H_parity_Block.block(offset_blocks, offset_blocks, M.rows(), M.cols()) = M;
        offset_blocks += M.rows();
    }
#endif
#ifdef showMatrix
        std::cout << H_parity_Block << "\n";
#endif
#ifdef saveMatrix
        saveComplexMatrixToFile(H_parity_Block, "HamiltonParityStates.txt", "parityStateAnsatz states Ansatz für N = " + std::to_string(N) +
                                                               "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), N);
#endif
#ifdef showEigenvalues
        std::cout << "eigenvalues:\n";
    for (std::complex<double> ev : *eiVals) {
        std::cout << ev << "\n";
    }
#endif
#ifdef saveEigenvalues
        saveComplexEiVals("EigenvaluesMParityStates.txt", "parityStateAnsatz states Ansatz für N = " + std::to_string(N) +
                                                       "\nJ1 = " + std::to_string(J1) + "\nJ2 = " + std::to_string(J2), *eiVals, N);
#endif
    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE) {
        const clock_t begin_time_PARITY = clock();

        std::cout << "\nparity states:..." << std::endl;

        auto *parityEiVals = new std::vector<std::complex<double>>;
        auto *matrixParityBlocks = new std::vector<Eigen::MatrixXcd>;

        getEiVals(J1, J2, parityEiVals, matrixParityBlocks, N, SIZE);

        auto time_PARITY = float(clock () - begin_time_PARITY) /  CLOCKS_PER_SEC;
        std::cout << "calculations done; this took: " << time_PARITY << " seconds\n";
        delete parityEiVals;
        delete matrixParityBlocks;
    }
}
