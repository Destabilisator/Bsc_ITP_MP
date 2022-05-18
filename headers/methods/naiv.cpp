#include "naiv.h"

/////////////////////////////// naiver Ansatz ///////////////////////////////

namespace naiv {
    void fillHamilton(double **hamilton, const double &J1,const double &J2, const int &N, const int &SIZE) {

        for (int s = 0; s <= SIZE - 1; s++) {
            for (int n = 0; n < N / 2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    hamilton[s][s] += 0.25 * J1;
                } else {
                    hamilton[s][s] -= 0.25 * J1;
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    hamilton[s][d] = 0.5 * J1;
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    hamilton[s][s] += 0.25 * J2;
                } else {
                    hamilton[s][s] -= 0.25 * J2;
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    hamilton[s][d] = 0.5 * J2;
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    hamilton[s][s] += 0.25 * J2;
                } else {
                    hamilton[s][s] -= 0.25 * J2;
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    hamilton[s][d] = 0.5 * J2;
                }
            }
        }

    }

    void getEiVals(const double &J1, const double &J2, std::vector<std::complex<double>> *HEiValList,
                   const int &N, const int &SIZE, Eigen::MatrixXcd &Matrix_U) {

        static auto **hamilton1 = new double*[SIZE];
        for (int i = 0; i < SIZE; i++) {
            hamilton1[i] = new double[SIZE];
            for (int j = 0; j < SIZE; j++) {
                hamilton1[i][j] = 0.0;
            }
        }
        fillHamilton(hamilton1, J1, J2, N, SIZE);
        Eigen::MatrixXd H1(SIZE, SIZE);
        for (int i = 0; i < SIZE; i++) {
            H1.row(i) = Eigen::VectorXd::Map(&hamilton1[i][0], SIZE);
        }
#ifdef showMatrix
        std::cout << H1 << std::endl;
#endif
#ifdef saveMatrix
        saveHamilton(hamilton1, "HamiltonNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), SIZE, N);
#endif
        //std::cout << "solving...\n";
        Eigen::EigenSolver<Eigen::MatrixXd> solver(H1);
        const Eigen::VectorXcd &H1EiVal = solver.eigenvalues();
        for (std::complex<double> ev : H1EiVal) {
            HEiValList->push_back(ev);
        }

        Matrix_U = solver.eigenvectors();

        HEiValList->shrink_to_fit();
#if defined(showEigenvalues) || defined(saveEigenvalues)
        HEiValList->shrink_to_fit();
        std::sort(HEiValList->begin(), HEiValList->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });
#endif
#ifdef showEigenvalues
        std::cout << "eigenvalues:\n";
        for (std::complex<double> ev : *HEiValList) {
            std::cout << ev << "\n";
        }
#endif
#ifdef saveEigenvalues
        saveComplexEiVals("EigenvaluesNaiv.txt", "naiver Ansatz für N = " + std::to_string(N), *HEiValList, N);
#endif

    }

    void start(const double &J1, const double &J2, const int &N, const int &SIZE, const double &START,
               const double &END, const int &COUNT, const int &cores) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "naiver Ansatz:..." << std::endl;
        auto *H_naiv_EiVals = new std::vector<std::complex<double>>;
        Eigen::MatrixXcd Matrix_U(SIZE, SIZE);
        getEiVals(J1, J2, H_naiv_EiVals, N, SIZE, Matrix_U);

        ///// susceptibility /////

        auto *susceptibility_naiv = new std::vector<std::tuple<double, double>>;

        Eigen::MatrixXd S2 = spinMatrix(N, SIZE);
        Eigen::MatrixXcd Matrix_U_inv_S2_U = Eigen::MatrixXcd::Zero(SIZE, SIZE);
        Matrix_U_inv_S2_U = Matrix_U.adjoint() * S2 * Matrix_U;

        for (int i = 0; i <= COUNT; i++) {
            double current = START + (END - START) * i / COUNT;
            //current_beta = 1 / current_beta;
            susceptibility_naiv->push_back({current,
                                            getSusceptibility(current, Matrix_U_inv_S2_U, *H_naiv_EiVals, N)});
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

        ///// save /////

        std::string filenameSusceptibility_X = "naiv_susceptibility.txt";
        std::string headerSusceptibility_X = "N: " + std::to_string(N) + "\n"
                                             + "T START: " + std::to_string(START) + "\n"
                                             + "T END: " + std::to_string(END) + "\n"
                                             + "data-points: " + std::to_string(COUNT) + "\n"
                                             + "calculation time with " + std::to_string(cores) + " threads: " + formatTime(elapsed_seconds);

        std::string headerWithJSusceptibility_X = "J1/J2 = " + std::to_string(J1/J2) +"\n" + headerSusceptibility_X;
        saveOutData(filenameSusceptibility_X, headerWithJSusceptibility_X, "J1/J2", "specific heat in J2", *susceptibility_naiv, N);
        std::cout << "\n";
    }

}
