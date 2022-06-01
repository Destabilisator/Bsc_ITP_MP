#include "multithreading.h"

/////////////////////////////// multi-threading ///////////////////////////////

namespace ED::multi {

    /////////////////////////////// Delta E, C(J) ///////////////////////////////
/*
    void get_DeltaE_CT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                             double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                             const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
        std::cout.flush();
        nextJMutex.unlock();

        while (pos <= COUNT) {

            std::vector<std::complex<double>> eiVals;
            std::vector<Eigen::MatrixXcd> matrixBlocks;

            momentumStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);

            // sort eigenvalues
            eiVals.shrink_to_fit();
            std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
                return std::real(c1) < std::real(c2);
            });

            // Delta E
            double E0 = std::real(eiVals.at(0));
            double E1 = std::real(eiVals.at(1));
            double specificHeat = getSpecificHeat(T, eiVals, N);

            // progressbar
            nextJMutex.lock();
            int prg = std::min({CURRENT, COUNT});
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = "
                        << START + (END - START) * prg / COUNT << " (" << prg << "/" << COUNT << ")     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outDataDeltaE->push_back({J, E1 - E0});
            outDataSpecificHeat_C->push_back({J, specificHeat});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

//            if (pos > COUNT) {
//                break;
//            } else {
//                J = START + (END - START) * pos / COUNT;
//            }

            J = START + (END - START) * pos / COUNT;

        }

    }

    void get_DeltaE_CT_const_parity(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                             double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                             const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
        std::cout.flush();
        nextJMutex.unlock();

        while (pos <= COUNT) {

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

            // progressbar
            nextJMutex.lock();
            int prg = std::min({CURRENT, COUNT});
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = "
                        << START + (END - START) * prg / COUNT << " (" << prg << "/" << COUNT << ")     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outDataDeltaE->push_back({J, E1 - E0});
            outDataSpecificHeat_C->push_back({J, specificHeat});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

//            if (pos > COUNT) {
//                break;
//            } else {
//                J = START + (END - START) * pos / COUNT;
//            }

            J = START + (END - START) * pos / COUNT;

        }

    }

    void get_DeltaE_CT_const_SI(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                    double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                    const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
        std::cout.flush();
        nextJMutex.unlock();

        while (pos <= COUNT) {

            std::vector<double> eiVals;
            std::vector<Eigen::MatrixXd> matrixBlocks;

            spinInversion::getEiVals(J, 1.0, &eiVals, &matrixBlocks, N, SIZE);

            // sort eigenvalues
            eiVals.shrink_to_fit();
            std::sort(eiVals.begin(), eiVals.end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
                return std::real(c1) < std::real(c2);
            });

            // Delta E
            double E0 = std::real(eiVals.at(0));
            double E1 = std::real(eiVals.at(1));
            double specificHeat = getSpecificHeat(T, eiVals, N);

            // progressbar
            nextJMutex.lock();
            int prg = std::min({CURRENT, COUNT});
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = "
                        << START + (END - START) * prg / COUNT << " (" << prg << "/" << COUNT << ")     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outDataDeltaE->push_back({J, E1 - E0});
            outDataSpecificHeat_C->push_back({J, specificHeat});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

//            if (pos > COUNT) {
//                break;
//            } else {
//                J = START + (END - START) * pos / COUNT;
//            }

            J = START + (END - START) * pos / COUNT;

        }

    }
*/

    void get_DeltaE_CT_const(double J, std::vector<std::tuple<double, double>> *outDataDeltaE,
                             double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                             const int &N, const int &SIZE) {

        std::vector<std::complex<double>> eiVals;
        std::vector<Eigen::MatrixXcd> matrixBlocks;

        momentumStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);

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

    void get_DeltaE_CT_const_SI(double J, std::vector<std::tuple<double, double>> *outDataDeltaE,
                                double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                                const int &N, const int &SIZE) {

        std::vector<double> eiVals;
        std::vector<Eigen::MatrixXd> matrixBlocks;

        spinInversion::getEiVals(J, 1.0, &eiVals, &matrixBlocks, N, SIZE);

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



    void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END,
                               int cores, const double &T, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "Delta E and C (T=const): calculating:...";

        std::vector<std::tuple<double, double>> outDataDeltaE;
        std::vector<std::tuple<double, double>> outDataSpecificHeat_C;

        std::cout << " (" << cores << ") cores";

        if (N%4 == 0) {
            if (N >= 12) {
                std::cout << ", spin inversion\n";
                int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataDeltaE, outDataSpecificHeat_C, T, N, SIZE, coutMutex, std::cout, curr)
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
                    std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << J << " (" << curr
                              << "/" << COUNT << ")     ";
                    std::cout.flush();
                    curr++;
                    coutMutex.unlock();

                    get_DeltaE_CT_const(J, &outDataDeltaE, T, &outDataSpecificHeat_C, N, SIZE);
                }
            } else {
                std::cout << ", parity states\n";
                int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataDeltaE, outDataSpecificHeat_C, T, N, SIZE, coutMutex, std::cout, curr)
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
                    std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << J << " (" << curr
                              << "/" << COUNT << ")     ";
                    std::cout.flush();
                    curr++;
                    coutMutex.unlock();

                    get_DeltaE_CT_const_parity(J, &outDataDeltaE, T, &outDataSpecificHeat_C, N, SIZE);
                }
            }

        } else {
            std::cout << ", momentum states\n";
            int curr = 0;
#pragma omp parallel for default(none) shared(START, END, COUNT, outDataDeltaE, outDataSpecificHeat_C, T, N, SIZE, coutMutex, std::cout, curr)
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
                std::cout << "] " << int((float) curr / (float) COUNT * 100) << "% J1/J2 = " << J << " (" << curr
                          << "/" << COUNT << ")     ";
                std::cout.flush();
                curr++;
                coutMutex.unlock();

                get_DeltaE_CT_const_SI(J, &outDataDeltaE, T, &outDataSpecificHeat_C, N, SIZE);
            }
        }



//        if (COUNT < cores) {
//            cores = COUNT;
//        }
//
//        cores = 1;
//
//        std::cout << " (" << cores << ") cores";
//
//        std::thread Threads[cores];
//
//        CURRENT = 0 + cores;
//
//        if (N%4 == 0) {
//            if (N >= 12) {
//                std::cout << ", spin inversion\n";
//                for (int i = 0; i < cores; i++) {
//                    Threads[i] = std::thread(get_DeltaE_CT_const_SI, START + (END - START) * i / COUNT, i + 1, &outDataDeltaE, T,
//                                             &outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
//                }
//            } else {
//                std::cout << ", parity states\n";
//                for (int i = 0; i < cores; i++) {
//                    Threads[i] = std::thread(get_DeltaE_CT_const_parity, START + (END - START) * i / COUNT, i + 1, &outDataDeltaE, T,
//                                             &outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
//                }
//            }
//
//        } else {
//            std::cout << ", momentum states\n";
//            for (int i = 0; i < cores; i++) {
//                Threads[i] = std::thread(get_DeltaE_CT_const, START + (END - START) * i / COUNT, i + 1, &outDataDeltaE, T,
//                                         &outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
//            }
//        }
//
//        for (int i = 0; i < cores; i++) {
//            Threads[i].join();
//        }

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
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: " +
                             std::to_string(elapsed_seconds.count()) + " seconds";

        std::string headerWithBeta = "T = " + std::to_string(T) + "\n" + header;

        saveOutData(filenameDeltaE, header, "J1/J2", "Delta E in J2", outDataDeltaE, N);
        saveOutData(filenameSpecificHeat_C, headerWithBeta, "J1/J2", "specific heat in J2", outDataSpecificHeat_C, N);
//        std::cout << "\n";

    }

    /////////////////////////////// Chi(J) ///////////////////////////////

    void get_XT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                      const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE, const double &T) {

        // progressbar init
        nextJMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
        std::cout.flush();
        nextJMutex.unlock();

        while (true) {

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
            // progressbar
            nextJMutex.lock();
            int prg = std::min({CURRENT, COUNT});
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = "
                        << START + (END - START) * prg / COUNT << " (" << prg << "/" << COUNT << ")     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outDataMagneticSusceptibility_X->push_back({J, susceptibility});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                break;
            } else {
                J = (double) START + (double) (END - START) * (double) pos / (double) COUNT;
            }
        }

    }

    void get_XT_const_momentum(double J, int pos, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                      const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE, const double &T) {

        // progressbar init
        nextJMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
        std::cout.flush();
        nextJMutex.unlock();

        while (true) {

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
            // progressbar
            nextJMutex.lock();
            int prg = std::min({CURRENT, COUNT});
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = "
                        << START + (END - START) * prg / COUNT << " (" << prg << "/" << COUNT << ")     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outDataMagneticSusceptibility_X->push_back({J, susceptibility});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                break;
            } else {
                J = (double) START + (double) (END - START) * (double) pos / (double) COUNT;
            }
        }

    }

    void start_XT_const(const int &COUNT, const double &START, const double &END,
                        int &cores, const double &T, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "X (T=const): calculating:...";

        std::vector<std::tuple<double, double>> outDataMagneticSusceptibility_X;

        if (COUNT < cores) {
            cores = COUNT;
        }
        std::cout << " (" << cores << ") cores";

        std::thread Threads[cores];

        CURRENT = 0 + cores;


        if (N >= 10) {
            std::cout << ", momentum states\n";
            for (int i = 0; i < cores; i++) {
                Threads[i] = std::thread(get_XT_const_momentum, START + (END - START) * i / COUNT, i + 1, &outDataMagneticSusceptibility_X, COUNT, START, END, N, SIZE, T);
            }
        } else {
            std::cout << ", magnetization blocks\n";
            for (int i = 0; i < cores; i++) {
                Threads[i] = std::thread(get_XT_const, START + (END - START) * i / COUNT, i + 1, &outDataMagneticSusceptibility_X, COUNT, START, END, N, SIZE, T);
            }
        }

        for (int i = 0; i < cores; i++) {
            Threads[i].join();
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

        std::string headerWithBeta = "T = " + std::to_string(T) + "\n" + header;

        saveOutData(filenameMagneticSusceptibility_X, headerWithBeta, "J1/J2", "magnetic susceptibility in J2", outDataMagneticSusceptibility_X, N);
//        std::cout << "\n";

    }

    /////////////////////////////// spin gap ///////////////////////////////

    void get_SpinGap_momentum(double J, int pos, std::vector<std::tuple<double, double>> *outDataSpinGap,
                              const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
        std::cout.flush();
        nextJMutex.unlock();

        while (true) {

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

            // progressbar
            nextJMutex.lock();
            int prg = std::min({CURRENT, COUNT});
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = "
                        << START + (END - START) * prg / COUNT << " (" << prg << "/" << COUNT << ")     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outDataSpinGap->push_back({J, E1-E0});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                break;
            } else {
                J = (double) START + (double) (END - START) * (double) pos / (double) COUNT;
            }
        }

    }

    void get_SpinGap_momentum_with_index(double J, int pos, std::vector<std::tuple<double, double, int, int, int, int>> *outDataSpinGap,
                              const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
        std::cout.flush();
        nextJMutex.unlock();

        while (true) {

//            std::vector<std::tuple<std::complex<double>, int, int>> data0;
//            std::vector<std::tuple<std::complex<double>, int, int>> data1;
//            momentumStates::getEiValsMagBlock_with_index(J, 1.0, data0, N, SIZE, N/2);
//            momentumStates::getEiValsMagBlock_with_index(J, 1.0, data1, N, SIZE, N/2 + 1);
//
//            // sort eigenvalues
//            data0.shrink_to_fit();
//            std::sort(data0.begin(), data0.end(), [](const std::tuple<std::complex<double>, int, int> &t1, const std::tuple<std::complex<double>, int, int> &t2) {
//                return std::real(std::get<0>(t1)) < std::real(std::get<0>(t2));
//            });
//
//            data1.shrink_to_fit();
//            std::sort(data1.begin(), data1.end(), [](const std::tuple<std::complex<double>, int, int> &t1, const std::tuple<std::complex<double>, int, int> &t2) {
//                return std::real(std::get<0>(t1)) < std::real(std::get<0>(t2));
//            });

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
//                    std::cout << E << " " << k << " " << S << std::endl;
                    // if nan, set it to first val
                    if (S == 0 && std::isnan(E0)) {E0 = E; k0 = k, S0 = S;}
                    else if ((S == 1 || S == -1) && std::isnan(E1)) {E1 = E; k1 = k, S1 = S;}
                    // if not nan check if E is lower
                    if (S == 0 && E < E0) {E0 = E; k0 = k, S0 = S;}
                    else if ((S == 1 || S == -1) && E < E1) {E1 = E; k1 = k, S1 = S;}
//                    else {continue;}
//                    std::cout << E0 << " " << k0 << " " << S0 << std::endl;
//                    std::cout << E1 << " " << k1 << " " << S1 << std::endl << std::endl;
                }
            }
//            std::cout << std::endl;
//            std::cout << E0 << " " << k0 << " " << S0 << std::endl;
//            std::cout << E1 << " " << k1 << " " << S1 << std::endl << std::endl;
//
//            double E0 = std::real(std::get<0>(data0.at(0)));
//            double E1 = std::real(std::get<0>(data1.at(0)));
//            int k0 = std::get<1>(data0.at(0));
//            int k1 = std::get<1>(data1.at(0));
//            int S0 = std::get<2>(data0.at(0));
//            int S1 = std::get<2>(data1.at(0));

            // progressbar
            nextJMutex.lock();
            int prg = std::min({CURRENT, COUNT});
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = "
                        << START + (END - START) * prg / COUNT << " (" << prg << "/" << COUNT << ")     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outDataSpinGap->push_back({J, E1-E0, k0, k1, S0, S1});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                break;
            } else {
                J = (double) START + (double) (END - START) * (double) pos / (double) COUNT;
            }
        }

    }

    void get_SpinGap_SI(double J, int pos, std::vector<std::tuple<double, double>> *outDataSpinGap,
                              const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        // progressbar init
        nextJMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
        std::cout.flush();
        nextJMutex.unlock();

        while (true) {

            std::vector<double> eiVals0;
            std::vector<double> eiVals1;
            spinInversion::getEiValsMagBlock(J, 1.0, &eiVals0, N, SIZE, N/2);
            spinInversion::getEiValsMagBlock(J, 1.0, &eiVals1, N, SIZE, N/2 + 1);
            spinInversion::getEiValsMagBlock(J, 1.0, &eiVals1, N, SIZE, N/2 - 1);

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

            // progressbar
            nextJMutex.lock();
            int prg = std::min({CURRENT, COUNT});
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = "
                        << START + (END - START) * prg / COUNT << " (" << prg << "/" << COUNT << ")     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outDataSpinGap->push_back({J, E1-E0});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                break;
            } else {
                J = (double) START + (double) (END - START) * (double) pos / (double) COUNT;
            }
        }

    }

    void start_SpinGap(const int &COUNT, const double &START, const double &END,
                       int &cores, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "spin gap: calculating:...";

        std::vector<std::tuple<double, double>> outDataSpinGap;

        if (COUNT < cores) {
            cores = COUNT;
        }
        std::cout << " (" << cores << ") cores";

        std::thread Threads[cores];

        CURRENT = 0 + cores;

        if (N%4 == 0) {
            std::cout << ", spin inversion\n";
            for (int i = 0; i < cores; i++) {
                Threads[i] = std::thread(get_SpinGap_SI, START + (END - START) * i / COUNT, i + 1, &outDataSpinGap, COUNT, START, END, N, SIZE);
            }

        } else {
            std::cout << ", momentum states\n";
            for (int i = 0; i < cores; i++) {
                Threads[i] = std::thread(get_SpinGap_momentum, START + (END - START) * i / COUNT, i + 1, &outDataSpinGap, COUNT, START, END, N, SIZE);
            }

        }

        for (int i = 0; i < cores; i++) {
            Threads[i].join();
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

        if (COUNT < cores) {
            cores = COUNT;
        }
        std::cout << " (" << cores << ") cores";

        std::thread Threads[cores];

        CURRENT = 0 + cores;

        std::cout << ", momentum states\n";
        for (int i = 0; i < cores; i++) {
            Threads[i] = std::thread(get_SpinGap_momentum_with_index, START + (END - START) * i / COUNT, i + 1, &outDataSpinGap, COUNT, START, END, N, SIZE);
        }


        for (int i = 0; i < cores; i++) {
            Threads[i].join();
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        // sort data-points
//        std::cout << "sorting final data" << std::endl;
        std::sort(outDataSpinGap.begin(), outDataSpinGap.end(), [](const std::tuple<double, double, int, int, int, int> &a, const std::tuple<double, double, int, int, int, int> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });
//        std::cout << "sorted final data" << std::endl;

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

}