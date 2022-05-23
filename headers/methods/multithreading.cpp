#include "multithreading.h"

/////////////////////////////// multi-threading ///////////////////////////////

namespace multi {

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

        while (true) {

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

            if (pos > COUNT) {
                break;
            } else {
                J = START + (END - START) * pos / COUNT;
            }

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

        while (true) {

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

            if (pos > COUNT) {
                break;
            } else {
                J = START + (END - START) * pos / COUNT;
            }

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

        while (true) {

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

            if (pos > COUNT) {
                break;
            } else {
                J = START + (END - START) * pos / COUNT;
            }

        }

    }

    void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END,
                               int &cores, const double &T, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "Delta E and C (T=const): calculating:...";

        std::vector<std::tuple<double, double>> outDataDeltaE;
        std::vector<std::tuple<double, double>> outDataSpecificHeat_C;

        if (COUNT < cores) {
            cores = COUNT;
        }
        std::cout << " (" << cores << ") cores";

        std::thread Threads[cores];

        CURRENT = 0 + cores;

        if (N%4 == 0) {
            if (N >= 12) {
                std::cout << ", spin inversion\n";
                for (int i = 0; i < cores; i++) {
                    Threads[i] = std::thread(get_DeltaE_CT_const_SI, START + (END - START) * i / COUNT, i + 1, &outDataDeltaE, T,
                                             &outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
                }
            } else {
                std::cout << ", parity states\n";
                for (int i = 0; i < cores; i++) {
                    Threads[i] = std::thread(get_DeltaE_CT_const_parity, START + (END - START) * i / COUNT, i + 1, &outDataDeltaE, T,
                                             &outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
                }
            }

        } else {
            std::cout << ", momentum states\n";
            for (int i = 0; i < cores; i++) {
                Threads[i] = std::thread(get_DeltaE_CT_const, START + (END - START) * i / COUNT, i + 1, &outDataDeltaE, T,
                                         &outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
            }
        }

        for (int i = 0; i < cores; i++) {
            Threads[i].join();
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
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: " +
                             std::to_string(elapsed_seconds.count()) + " seconds";

        std::string headerWithBeta = "T = " + std::to_string(T) + "\n" + header;

        saveOutData(filenameDeltaE, header, "J1/J2", "Delta E in J2", outDataDeltaE, N);
        saveOutData(filenameSpecificHeat_C, headerWithBeta, "J1/J2", "specific heat in J2", outDataSpecificHeat_C, N);

        std::cout << "\n";

    }

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

        std::cout << "\n";

    }

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

    void get_SpinGap_momentum_with_k(double J, int pos, std::vector<std::tuple<double, double, int, int>> *outDataSpinGap,
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

            std::vector<std::tuple<std::complex<double>, int>> data0;
            std::vector<std::tuple<std::complex<double>, int>> data1;
            momentumStates::getEiValsMagBlock_with_k(J, 1.0, data0, N, SIZE, N/2);
            momentumStates::getEiValsMagBlock_with_k(J, 1.0, data1, N, SIZE, N/2 + 1);

            // sort eigenvalues
            data0.shrink_to_fit();
            std::sort(data0.begin(), data0.end(), [](const std::tuple<std::complex<double>, int> &t1, const std::tuple<std::complex<double>, int> &t2) {
                return std::real(std::get<0>(t1)) < std::real(std::get<0>(t2));
            });

            data1.shrink_to_fit();
            std::sort(data1.begin(), data1.end(), [](const std::tuple<std::complex<double>, int> &t1, const std::tuple<std::complex<double>, int> &t2) {
                return std::real(std::get<0>(t1)) < std::real(std::get<0>(t2));
            });

            // Delta E
            double E0 = std::real(std::get<0>(data0.at(0)));
            double E1 = std::real(std::get<0>(data1.at(0)));

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
            outDataSpinGap->push_back({J, E1-E0, std::get<1>(data0.at(0)), std::get<1>(data1.at(0))});
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

        std::string filenameMagneticSusceptibility_X = "data_spin_gap.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: "
                             + std::to_string(elapsed_seconds.count()) + " seconds";

        saveOutData(filenameMagneticSusceptibility_X, header, "J1/J2", "spin gap in J2", outDataSpinGap, N);

        std::cout << "\n";

    }

    void start_SpinGap_with_k(const int &COUNT, const double &START, const double &END,
                       int &cores, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "spin gap: calculating:...";

        std::vector<std::tuple<double, double, int, int>> outDataSpinGap;

        if (COUNT < cores) {
            cores = COUNT;
        }
        std::cout << " (" << cores << ") cores";

        std::thread Threads[cores];

        CURRENT = 0 + cores;

        std::cout << ", momentum states\n";
        for (int i = 0; i < cores; i++) {
            Threads[i] = std::thread(get_SpinGap_momentum_with_k, START + (END - START) * i / COUNT, i + 1, &outDataSpinGap, COUNT, START, END, N, SIZE);
        }


        for (int i = 0; i < cores; i++) {
            Threads[i].join();
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        // sort data-points
        std::sort(outDataSpinGap.begin(), outDataSpinGap.end(), [](const std::tuple<double, double, int, int> &a, const std::tuple<double, double, int, int> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::string filenameMagneticSusceptibility_X = "data_spin_gap.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: "
                             + std::to_string(elapsed_seconds.count()) + " seconds";

        saveOutData(filenameMagneticSusceptibility_X, header, "J1/J2", "spin gap in J2\tk1\tk2", outDataSpinGap, N);

        std::cout << "\n";

    }

}