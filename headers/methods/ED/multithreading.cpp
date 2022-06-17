#include "multithreading.h"

/////////////////////////////// multi-threading ///////////////////////////////

namespace ED::multi {

    /////////////////////////////// Delta E, C(J) ///////////////////////////////

    void get_DeltaE_CT_const_momentum(double J, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                      double temp, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                      const double &h, const int &N, const int &SIZE) {

        std::vector<std::complex<double>> eiVals;
        std::vector<Eigen::MatrixXcd> matrixBlocks;

        momentumStates::getEiVals(J, 1.0, h, eiVals, matrixBlocks, N, SIZE);

        // sort eigenvalues
        eiVals.shrink_to_fit();
        std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });

        // Delta E
        double E0 = std::real(eiVals.at(0));
        double E1 = std::real(eiVals.at(1));
        double specificHeat = getSpecificHeat(temp, eiVals, N);


        // write data
        nextJMutex.lock();
        outDataDeltaE->push_back({J, E1 - E0});
        outDataSpecificHeat_C->push_back({J, specificHeat});
//        CURRENT++;
        nextJMutex.unlock();

    }
/*
    void get_DeltaE_CT_const_parity(double J, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                    double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                    const int &N, const int &SIZE) {

        std::vector<double> eiVals;
        std::vector<Eigen::MatrixXd> matrixBlocks;

        parityStates::getEiVals(J, 1.0, &eiVals, &matrixBlocks, N, SIZE);

        // sort eigenvalues
        eiVals.shrink_to_fit();
        std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });

        // Delta E
        double E0 = std::real(eiVals.at(0));
        double E1 = std::real(eiVals.at(1));
        double specificHeat = getSpecificHeat(T, eiVals, N);

        // write data
        nextJMutex.lock();
        outDataDeltaE->push_back({J, E1 - E0});
        outDataSpecificHeat_C->push_back({J, specificHeat});
//        CURRENT++;
        nextJMutex.unlock();

    }
*/
    void get_DeltaE_CT_const_SI(double J, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                double temp, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                const double &h, const int &N, const int &SIZE) {

        std::vector<double> eiVals;
        std::vector<Eigen::MatrixXd> matrixBlocks;

        spinInversion::getEiVals(J, 1.0, h, &eiVals, &matrixBlocks, N, SIZE);

        // sort eigenvalues
        eiVals.shrink_to_fit();
        std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });

        // Delta E
        double E0 = std::real(eiVals.at(0));
        double E1 = std::real(eiVals.at(1));
        double specificHeat = getSpecificHeat(temp, eiVals, N);

        // write data
        nextJMutex.lock();
        outDataDeltaE->push_back({J, E1 - E0});
        outDataSpecificHeat_C->push_back({J, specificHeat});
//        CURRENT++;
        nextJMutex.unlock();

    }



    void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END, const double &h,
                               int cores, const double &temp, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "Delta E and C (T=const): calculating:...";

        std::vector<std::tuple<double, double>> outDataDeltaE;
        std::vector<std::tuple<double, double>> outDataSpecificHeat_C;

        std::cout << " (" << cores << ") cores";

        if (N%4 == 0) {
//            if (N >= 12) {
                std::cout << ", spin inversion\n";
                int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataDeltaE, outDataSpecificHeat_C, temp, N, SIZE, coutMutex, std::cout, curr, h)
                for (int i = 0; i <= COUNT; i++) {
                    double J = START + (END - START) * i / COUNT;
                    coutMutex.lock();
                    int p = (int) ((float) curr / (float) (COUNT) * (float) PROGRESSBAR_SEGMENTS);
                    std::cout << "\r[";
                    for (int _ = 0; _ < p; _++) {
                        std::cout << "#";
                    }
                    for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                        std::cout << ".";
                    }
                    std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << START + (END - START) * curr / COUNT << " (" << curr
                              << "/" << COUNT << ")     ";
                    std::cout.flush();
                    curr++;
                    coutMutex.unlock();

                    get_DeltaE_CT_const_SI(J, &outDataDeltaE, temp, &outDataSpecificHeat_C, h, N, SIZE);
                }
//            } else {
//                std::cout << ", parity states\n";
//                int curr = 0;
//#pragma omp parallel for default(none) shared(START, END, COUNT, outDataDeltaE, outDataSpecificHeat_C, T, N, SIZE, coutMutex, std::cout, curr)
//                for (int i = 0; i <= COUNT; i++) {
//                    double J = START + (END - START) * i / COUNT;
//                    coutMutex.lock();
//                    int p = (int) ((float) curr / (float) (COUNT) * (float) PROGRESSBAR_SEGMENTS);
//                    std::cout << "\r[";
//                    for (int _ = 0; _ < p; _++) {
//                        std::cout << "#";
//                    }
//                    for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
//                        std::cout << ".";
//                    }
//                    std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << START + (END - START) * curr / COUNT << " (" << curr
//                              << "/" << COUNT << ")     ";
//                    std::cout.flush();
//                    curr++;
//                    coutMutex.unlock();
//
//                    get_DeltaE_CT_const_parity(J, &outDataDeltaE, T, &outDataSpecificHeat_C, N, SIZE);
//                }
//            }
        } else {
            std::cout << ", momentum states\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataDeltaE, outDataSpecificHeat_C, temp, N, SIZE, coutMutex, std::cout, curr, h)
            for (int i = 0; i <= COUNT; i++) {
                double J = START + (END - START) * i / COUNT;
                coutMutex.lock();
                int p = (int) ((float) curr / (float) (COUNT) * (float) PROGRESSBAR_SEGMENTS);
                std::cout << "\r[";
                for (int _ = 0; _ < p; _++) {
                    std::cout << "#";
                }
                for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                    std::cout << ".";
                }
                std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << START + (END - START) * curr / COUNT << " (" << curr
                          << "/" << COUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_DeltaE_CT_const_momentum(J, &outDataDeltaE, temp, &outDataSpecificHeat_C, h, N, SIZE);
            }
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        // sort data-points
        std::sort(outDataDeltaE.begin(), outDataDeltaE.end(), [](
                const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });
        std::sort(outDataSpecificHeat_C.begin(), outDataSpecificHeat_C.end(), [](
                const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::string filenameDeltaE = "data_delta_E.txt";
        std::string filenameSpecificHeat_C = "data_specific_heat_T_const.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "h: " + std::to_string(h) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: " +
                             std::to_string(elapsed_seconds.count()) + " seconds";


        std::string headerWithBeta = "beta: " + std::to_string(temp) + "\n" + header;

        saveOutData(filenameDeltaE, header, "J1/J2", "Delta E in J2", outDataDeltaE, N);
        saveOutData(filenameSpecificHeat_C, headerWithBeta, "J1/J2", "specific heat in J2", outDataSpecificHeat_C, N);
//        std::cout << "\n";

    }

    /////////////////////////////// Chi(J) ///////////////////////////////

    void get_XT_const_magnetization(double J, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                                    const int &N, const int &SIZE, const double &T) {

        std::vector<int> states;
        std::vector<std::complex<double>> eiVals;

        fillStates(&states, N/2, N, SIZE);
        const int statesCount = (int) states.size();

        Eigen::MatrixXcd matrixBlockU;
        magnetizationBlocks::getEiValsZeroBlock(J, 1.0, &eiVals, matrixBlockU, states, N);

        Eigen::MatrixXd S2 = spinMatrix(N, states);
        Eigen::MatrixXcd Matrix_U_inv_S2_U = Eigen::MatrixXcd::Zero(statesCount, statesCount);
        Matrix_U_inv_S2_U = matrixBlockU.adjoint() * S2 * matrixBlockU;

        double susceptibility = getSusceptibilityDegeneracy(T, Matrix_U_inv_S2_U, eiVals, N);

        // write data
        nextJMutex.lock();
        outDataMagneticSusceptibility_X->push_back({J, susceptibility});
//        CURRENT++;
        nextJMutex.unlock();


    }

    void get_XT_const_momentum(double J, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                               const int &N, const int &SIZE, const double &T){

        std::vector<std::vector<std::complex<double>>> eiVals;
        std::vector<Eigen::MatrixXcd> matrixBlockU;
        std::vector<Eigen::MatrixXcd> matrixBlockS2;
        momentumStates::getEiValsZeroBlock(J, 1.0, eiVals, matrixBlockU, matrixBlockS2, N, SIZE);

        std::vector<Eigen::MatrixXcd> Blocks_U_inv_S2_U;
        for(int i = 0; i < matrixBlockU.size(); i++) {
            Eigen::MatrixXcd M = matrixBlockU.at(i).adjoint() * matrixBlockS2.at(i) * matrixBlockU.at(i);
            Blocks_U_inv_S2_U.push_back(M);
        }

        double susceptibility = getSusceptibilityDegeneracy(T, Blocks_U_inv_S2_U, eiVals, N);

        // write data
        nextJMutex.lock();
        outDataMagneticSusceptibility_X->push_back({J, susceptibility});
//            CURRENT++;
            nextJMutex.unlock();


    }

    void start_XT_const(const int &COUNT, const double &START, const double &END,
                        int &cores, const double &T, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "X (T=const): calculating:...";

        std::vector<std::tuple<double, double>> outDataMagneticSusceptibility_X;

        std::cout << " (" << cores << ") cores";

        if (N >= 10) {
            std::cout << ", momentum states\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataMagneticSusceptibility_X, T, N, SIZE, coutMutex, std::cout, curr)
            for (int i = 0; i <= COUNT; i++) {
                double J = START + (END - START) * i / COUNT;
                coutMutex.lock();
                int p = (int) ((float) curr / (float) (COUNT) * (float) PROGRESSBAR_SEGMENTS);
                std::cout << "\r[";
                for (int _ = 0; _ < p; _++) {
                    std::cout << "#";
                }
                for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                    std::cout << ".";
                }
                std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << START + (END - START) * curr / COUNT << " (" << curr
                          << "/" << COUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_XT_const_momentum(J, &outDataMagneticSusceptibility_X, N, SIZE, T);
            }
        } else {
            std::cout << ", magnetization blocks\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataMagneticSusceptibility_X, T, N, SIZE, coutMutex, std::cout, curr)
            for (int i = 0; i <= COUNT; i++) {
                double J = START + (END - START) * i / COUNT;
                coutMutex.lock();
                int p = (int) ((float) curr / (float) (COUNT) * (float) PROGRESSBAR_SEGMENTS);
                std::cout << "\r[";
                for (int _ = 0; _ < p; _++) {
                    std::cout << "#";
                }
                for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                    std::cout << ".";
                }
                std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << START + (END - START) * curr / COUNT << " (" << curr
                          << "/" << COUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_XT_const_magnetization(J, &outDataMagneticSusceptibility_X, N, SIZE, T);
            }
        }


        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        // sort data-points
        std::sort(outDataMagneticSusceptibility_X.begin(), outDataMagneticSusceptibility_X.end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::string filenameMagneticSusceptibility_X = "data_susceptibility_T_const.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: "
                             + std::to_string(elapsed_seconds.count()) + " seconds";

        std::string headerWithBeta = "T: " + std::to_string(T) + "\n" + header;

        saveOutData(filenameMagneticSusceptibility_X, headerWithBeta, "J1/J2", "magnetic susceptibility in J2", outDataMagneticSusceptibility_X, N);
//        std::cout << "\n";

    }

    /////////////////////////////// spin gap ///////////////////////////////

    void get_SpinGap_momentum(std::vector<std::tuple<double, double>> *outDataSpinGap,
                              const double &J, const int &N, const int &SIZE) {

        std::vector<std::complex<double>> data0;
        std::vector<std::complex<double>> data1;
        momentumStates::getEiValsMagBlock(J, 1.0, data0, N, SIZE, N/2);
        momentumStates::getEiValsMagBlock(J, 1.0, data1, N, SIZE, N/2 + 1);
        momentumStates::getEiValsMagBlock(J, 1.0, data1, N, SIZE, N/2 - 1);

        // sort eigenvalues
        data0.shrink_to_fit();
        std::sort(data0.begin(), data0.end(), [](const std::complex<double> &c1, const std::complex<double>&c2) {
            return std::real(c1) < std::real(c2);
        });

        data1.shrink_to_fit();
        std::sort(data1.begin(), data1.end(), [](const std::complex<double> &c1, const std::complex<double>&c2) {
            return std::real(c1) < std::real(c2);
        });

        // Delta E
        double E0 = std::real(data0.at(0));
        double E1 = std::real(data1.at(0));

        // write data
        nextJMutex.lock();
        outDataSpinGap->push_back({J, E1-E0});
//        CURRENT++;
        nextJMutex.unlock();

    }

    void get_SpinGap_momentum_with_index(std::vector<std::tuple<double, double, int, int, int, int>> *outDataSpinGap,
                                         const double &J, const int &N, const int &SIZE) {

        std::vector<std::vector<std::tuple<std::complex<double>, int, int>>> data(N+1);
        for (int m = 0; m <= N; m++) {
            std::vector<std::tuple<std::complex<double>, int, int>> d;
            momentumStates::getEiValsMagBlock_with_index(J, 1.0, d, N, SIZE, m);
            d.shrink_to_fit();
            std::sort(d.begin(), d.end(), [](const std::tuple<std::complex<double>, int, int> &t1, const std::tuple<std::complex<double>, int, int> &t2) {
                return std::real(std::get<0>(t1)) < std::real(std::get<0>(t2));
            });
            data.push_back(d);
        }

        double E0 = NAN;
        double E1 = NAN;
        int k0, k1, S0, S1;
        for (const std::vector<std::tuple<std::complex<double>, int, int>>& dat : data) {
            for (std::tuple<std::complex<double>, int, int> d : dat) {
                double E = std::real(std::get<0>(d));
                int k = std::get<1>(d);
                int S = std::get<2>(d);
                // if nan, set it to first val
                if (S == 0 && std::isnan(E0)) {E0 = E; k0 = k, S0 = S;}
                else if ((S == 1 || S == -1) && std::isnan(E1)) {E1 = E; k1 = k, S1 = S;}
                // if not nan check if E is lower
                if (S == 0 && E < E0) {E0 = E; k0 = k, S0 = S;}
                else if ((S == 1 || S == -1) && E < E1) {E1 = E; k1 = k, S1 = S;}
            }
        }
        // write data
        nextJMutex.lock();
        outDataSpinGap->push_back({J, E1-E0, k0, k1, S0, S1});
//        CURRENT++;
        nextJMutex.unlock();

    }

    void get_SpinGap_SI(std::vector<std::tuple<double, double>> *outDataSpinGap,
                        const double &J, const int &N, const int &SIZE) {

        std::vector<double> eiVals0;
        std::vector<double> eiVals1;
        spinInversion::getEiValsMagBlock(J, 1.0, 0.0, &eiVals0, N, SIZE, N/2);
        spinInversion::getEiValsMagBlock(J, 1.0, 0.0, &eiVals1, N, SIZE, N/2 + 1);
        spinInversion::getEiValsMagBlock(J, 1.0, 0.0, &eiVals1, N, SIZE, N/2 - 1);

        // sort eigenvalues
        eiVals0.shrink_to_fit();
        std::sort(eiVals0.begin(), eiVals0.end(), [](const double &d1, const double &d2) {
            return d1 < d2;
        });

        eiVals1.shrink_to_fit();
        std::sort(eiVals1.begin(), eiVals1.end(), [](const double &d1, const double &d2) {
            return d1 < d2;
        });

        // Delta E
        double E0 = eiVals0.at(0);
        double E1 = eiVals1.at(0);

        // write data
        nextJMutex.lock();
        outDataSpinGap->push_back({J, E1-E0});
//        CURRENT++;
        nextJMutex.unlock();

    }

    void start_SpinGap(const int &COUNT, const double &START, const double &END,
                       int &cores, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "spin gap: calculating:...";

        std::vector<std::tuple<double, double>> outDataSpinGap;

        std::cout << " (" << cores << ") cores";

        if (N%4 == 0) {
            std::cout << ", spin inversion\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataSpinGap, N, SIZE, coutMutex, std::cout, curr)
            for (int i = 0; i <= COUNT; i++) {
                double J = START + (END - START) * i / COUNT;
                coutMutex.lock();
                int p = (int) ((float) curr / (float) (COUNT) * (float) PROGRESSBAR_SEGMENTS);
                std::cout << "\r[";
                for (int _ = 0; _ < p; _++) {
                    std::cout << "#";
                }
                for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                    std::cout << ".";
                }
                std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << START + (END - START) * curr / COUNT << " (" << curr
                          << "/" << COUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_SpinGap_SI(&outDataSpinGap, J, N, SIZE);
            }
        } else {
            std::cout << ", momentum states\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataSpinGap, N, SIZE, coutMutex, std::cout, curr)
            for (int i = 0; i <= COUNT; i++) {
                double J = START + (END - START) * i / COUNT;
                coutMutex.lock();
                int p = (int) ((float) curr / (float) (COUNT) * (float) PROGRESSBAR_SEGMENTS);
                std::cout << "\r[";
                for (int _ = 0; _ < p; _++) {
                    std::cout << "#";
                }
                for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                    std::cout << ".";
                }
                std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << START + (END - START) * curr / COUNT << " (" << curr
                          << "/" << COUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_SpinGap_momentum(&outDataSpinGap, J, N, SIZE);
            }
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        // sort data-points
        std::sort(outDataSpinGap.begin(), outDataSpinGap.end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::string filenameMagneticSusceptibility_X = "data_spin_gap.txt"; // data_spin_gap
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: "
                             + std::to_string(elapsed_seconds.count()) + " seconds";

        saveOutData(filenameMagneticSusceptibility_X, header, "J1/J2", "spin gap in J2", outDataSpinGap, N);
//        std::cout << "\n";

    }

    void start_SpinGap_with_index(const int &COUNT, const double &START, const double &END,
                       int &cores, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "spin gap: calculating:...";

        std::vector<std::tuple<double, double, int, int, int, int>> outDataSpinGap;

        std::cout << " (" << cores << ") cores";

        std::cout << ", momentum states\n";
        int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataSpinGap, N, SIZE, coutMutex, std::cout, curr)
        for (int i = 0; i <= COUNT; i++) {
            double J = START + (END - START) * i / COUNT;
            coutMutex.lock();
            int p = (int) ((float) curr / (float) (COUNT) * (float) PROGRESSBAR_SEGMENTS);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            }
            for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            }
            std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << START + (END - START) * curr / COUNT << " (" << curr
                      << "/" << COUNT << ")     ";
            std::cout.flush();
            curr++;
            coutMutex.unlock();

            get_SpinGap_momentum_with_index(&outDataSpinGap, J, N, SIZE);
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        // sort data-points
        std::sort(outDataSpinGap.begin(), outDataSpinGap.end(), [](const std::tuple<double, double, int, int, int, int> &a, const std::tuple<double, double, int, int, int, int> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });
        std::string filenameMagneticSusceptibility_X = "data_spin_gap_with_index_V2.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: "
                             + std::to_string(elapsed_seconds.count()) + " seconds";

        saveOutData(filenameMagneticSusceptibility_X, header, "J1/J2", "spin gap in J2\tk0\tk1\tS0\tS1", outDataSpinGap, N);
//        std::cout << "\n";

    }

    void startSusceptibilityMultiJ(const double &J_START, const double &J_END, const int &J_COUNT,
                                   const double &BETA_START, const double &BETA_END, const double &BETA_COUNT,
                                   const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "susceptibility (momentum states): calculating..." << std::endl;

        ///// gather x-data /////
        std::vector<double> beta_Data;
        for (int i = 0; i <= BETA_COUNT; i++) {
            double current = BETA_START + (BETA_END - BETA_START) * i / BETA_COUNT;
            beta_Data.emplace_back(current);
        } beta_Data.shrink_to_fit();

        // init progressbar
        int prgbar_segm = 50;
        int curr = 0;
        coutMutex.lock();
        int _p = (int) ( (float) curr / (float) J_COUNT * (float) prgbar_segm);
        std::cout << "\r[";
        for (int _ = 0; _ < _p; _++) {
            std::cout << "#";
        } for (int _ = _p; _ < prgbar_segm; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) curr / (float) J_COUNT * 100.0 ) << "% (" << curr << "/" << J_COUNT << "), J1/J2 = " << J_START + (J_END - J_START) * curr / J_COUNT << "     ";
        std::cout.flush();
        curr++;
        coutMutex.unlock();

        ///// susceptibility /////

#pragma omp parallel for default(none) shared(J_COUNT, J_START, J_END, N, SIZE, coutMutex, BETA_START, BETA_END, BETA_COUNT, curr, prgbar_segm, std::cout, beta_Data)
        for (int J_pos = 0; J_pos < J_COUNT; J_pos++) {

            double J = J_START + (J_END - J_START) * J_pos / J_COUNT;

            std::vector<std::vector<std::complex<double>>> eiVals;
            std::vector<Eigen::MatrixXcd> matrixBlockU;
            std::vector<Eigen::MatrixXcd> matrixBlockS2;
            momentumStates::getEiValsZeroBlock(J, 1.0, eiVals, matrixBlockU, matrixBlockS2, N, SIZE);

            std::vector<std::tuple<double, double>> susceptibility_magnetization;

            std::vector<Eigen::MatrixXcd> Blocks_U_inv_S2_U;
            for (int i = 0; i < matrixBlockU.size(); i++) {
                Eigen::MatrixXcd M = matrixBlockU.at(i).adjoint() * matrixBlockS2.at(i) * matrixBlockU.at(i);
                Blocks_U_inv_S2_U.push_back(M);
            }

            for (double beta : beta_Data) {
                susceptibility_magnetization.emplace_back(beta, getSusceptibilityDegeneracy(beta, Blocks_U_inv_S2_U, eiVals, N));
            }

            ///// save /////

            saveOutDataSilent("/spin_gap_data/1/X_J" + std::to_string(J) + "ED.txt",
                              "\n", "J1/J2", "specific heat in J2", susceptibility_magnetization, N);

            // progressbar
            coutMutex.lock();
            int p = (int) ( (float) curr / (float) J_COUNT * (float) prgbar_segm);
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < prgbar_segm; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) curr / (float) J_COUNT * 100.0 ) << "% (" << curr << "/" << J_COUNT << "), J1/J2 = " << J_START + (J_END - J_START) * curr / J_COUNT << "     ";
            std::cout.flush();
            curr++;
            coutMutex.unlock();

        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

    }

//    /////////////////////////////// excitation energies ///////////////////////////////
//
//    void startSpecificHeatMultiJ(const double &J_START, const double &J_END, const int &J_COUNT,
//                                   const double &BETA_START, const double &BETA_END, const double &BETA_COUNT,
//                                   const int &N, const int &SIZE)

}