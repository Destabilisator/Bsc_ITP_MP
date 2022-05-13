#include "multithreading.h"

/////////////////////////////// multi-threading ///////////////////////////////

namespace multi {
    void get_DeltaE_CT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                             double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                             const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        while (true) {

            int p = (int) ( (float) pos / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) pos / (float) COUNT * 100.0 ) << "% J1/J2 = " << J << " (" << pos << "/" << COUNT << ")     ";
            coutMutex.unlock();

            auto *eiVals = new std::vector<std::complex<double>>;
            auto *matrixBlocks = new std::vector<Eigen::MatrixXcd>;
//            auto *eiVals = new std::vector<double>;
//            auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;

//            naiv::getEiVals(J, 1.0, eiVals, N, SIZE);
//            magnetizationBlocks::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);
            momentumStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);
//            parityStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);

            // sort eigenvalues
            eiVals->shrink_to_fit();
            std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
                return std::real(c1) < std::real(c2);
            });

            // Delta E
            double E0 = std::real(eiVals->at(0));
            double E1 = std::real(eiVals->at(1));

            // write data
            nextJMutex.lock();
            outDataDeltaE->push_back({J, E1 - E0});
            outDataSpecificHeat_C->push_back({J, getSpecificHeat(T, *eiVals, N)});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                delete eiVals;
                delete matrixBlocks;
                break;
            } else {
                J = START + (END - START) * pos / COUNT;
            }
            // clean up
            eiVals->clear();
            matrixBlocks->clear();

        }

    }

    void get_DeltaE_CT_const_parity(double J, int pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                             double T, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                             const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

        while (true) {

            int p = (int) ( (float) pos / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) pos / (float) COUNT * 100.0 ) << "% J1/J2 = " << J << " (" << pos << "/" << COUNT << ")     ";
            coutMutex.unlock();

//            auto *eiVals = new std::vector<std::complex<double>>;
//            auto *matrixBlocks = new std::vector<Eigen::MatrixXcd>;
            auto *eiVals = new std::vector<double>;
            auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;

//            naiv::getEiVals(J, 1.0, eiVals, N, SIZE);
//            magnetizationBlocks::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);
//            momentumStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);
            parityStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);

            // sort eigenvalues
            eiVals->shrink_to_fit();
            std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
                return std::real(c1) < std::real(c2);
            });

            // Delta E
            double E0 = std::real(eiVals->at(0));
            double E1 = std::real(eiVals->at(1));

            // write data
            nextJMutex.lock();
            outDataDeltaE->push_back({J, E1 - E0});
            outDataSpecificHeat_C->push_back({J, getSpecificHeat(T, *eiVals, N)});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                delete eiVals;
                delete matrixBlocks;
                break;
            } else {
                J = START + (END - START) * pos / COUNT;
            }
            // clean up
            eiVals->clear();
            matrixBlocks->clear();

        }

    }

    void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END, const unsigned int *cpu_cnt,
                               int *cores, const double &T, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "Delta E and C (T=const): calculating:..." << std::endl;

        auto *outDataDeltaE = new std::vector<std::tuple<double, double>>;
        auto *outDataSpecificHeat_C = new std::vector<std::tuple<double, double>>;

        if (COUNT < *cpu_cnt) {
            *cores = COUNT;
        }
        std::thread Threads[*cores];

        CURRENT += *cores;

        if (N%4 == 0) {
            for (int i = 0; i < *cores; i++) {
                Threads[i] = std::thread(get_DeltaE_CT_const_parity, START + (END - START) * i / COUNT, i + 1, outDataDeltaE, T,
                                         outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
            }
        } else {
            for (int i = 0; i < *cores; i++) {
                Threads[i] = std::thread(get_DeltaE_CT_const, START + (END - START) * i / COUNT, i + 1, outDataDeltaE, T,
                                         outDataSpecificHeat_C, COUNT, START, END, N, SIZE);
            }
        }



        for (int i = 0; i < *cores; i++) {
            Threads[i].join();
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << elapsed_seconds.count() << " seconds\n\n";

        // sort data-points
        std::sort(outDataDeltaE->begin(), outDataDeltaE->end(), [](
                const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });
        std::sort(outDataSpecificHeat_C->begin(), outDataSpecificHeat_C->end(), [](
                const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::string filenameDeltaE = "data_delta_E.txt";
        std::string filenameSpecificHeat_C = "data_specific_heat.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(*cores) + " threads: " +
                             std::to_string(elapsed_seconds.count()) + " seconds";

        std::string headerWithBeta = "T = " + std::to_string(T) + "\n" + header;

        saveOutData(filenameDeltaE, header, "J1/J2", "Delta E in J2", *outDataDeltaE, N);
        saveOutData(filenameSpecificHeat_C, headerWithBeta, "J1/J2", "specific heat in J2", *outDataSpecificHeat_C, N);

        delete outDataDeltaE;
        delete outDataSpecificHeat_C;

    }

    void get_XT_const(double J, int pos, std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X,
                      const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE, const double &T) {

        while (true) {

            int p = (int) ( (float) pos / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) pos / (float) COUNT * 100.0 ) << "% J1/J2 = " << J << " (" << pos << "/" << COUNT << ")     ";
            coutMutex.unlock();


            auto *states = new std::vector<int>;
            auto *eiVals = new std::vector<std::complex<double>>;
            Eigen::MatrixXcd matrixBlockU;
            magnetizationBlocks::getEiValsZeroBlock(J, 1.0, eiVals, matrixBlockU, states, N, SIZE);

            ///// susceptibility /////

            Eigen::MatrixXd S2 = spinMatrix(N, *states);
            Eigen::MatrixXcd Matrix_U_inv_S2_U = Eigen::MatrixXcd::Zero(SIZE, SIZE);
            Matrix_U_inv_S2_U = matrixBlockU.adjoint() * S2 * matrixBlockU;

            double susceptibility = getSusceptibilityDegeneracy(T, Matrix_U_inv_S2_U, *eiVals, N);

            // write data
            nextJMutex.lock();
            outDataMagneticSusceptibility_X->push_back({J, susceptibility});
            pos = CURRENT;
            CURRENT++;
            nextJMutex.unlock();

            if (pos > COUNT) {
                delete eiVals;
                delete states;
                break;
            } else {
                J = START + (END - START) * pos / COUNT;
            }
            // clean up
            eiVals->clear();
            states->clear();

        }

    }

    void start_XT_const(const int &COUNT, const double &START, const double &END, const unsigned int *cpu_cnt,
                        int *cores, const double &T, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "X (T=const): calculating:..." << std::endl;

        auto *outDataMagneticSusceptibility_X = new std::vector<std::tuple<double, double>>;

        if (COUNT < *cpu_cnt) {
            *cores = COUNT;
        }
        std::thread Threads[*cores];

        CURRENT = 1 + *cores;

        for (int i = 0; i < *cores; i++) {
            Threads[i] = std::thread(get_XT_const, START + (END - START) * i / COUNT, i + 1, outDataMagneticSusceptibility_X, COUNT, START, END, N, SIZE, T);
        }

        for (int i = 0; i < *cores; i++) {
            Threads[i].join();
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << elapsed_seconds.count() << " seconds\n\n";

        // sort data-points
        std::sort(outDataMagneticSusceptibility_X->begin(), outDataMagneticSusceptibility_X->end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::string filenameMagneticSusceptibility_X = "data_magnetic_susceptibility_T_const.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "J1/J2 START: " + std::to_string(START) + "\n"
                             + "J1/J2 END: " + std::to_string(END) + "\n"
                             + "data-points: " + std::to_string(COUNT) + "\n"
                             + "calculation time with " + std::to_string(*cores) + " threads: "
                             + std::to_string(elapsed_seconds.count()) + " seconds";

        std::string headerWithBeta = "T = " + std::to_string(T) + "\n" + header;

        saveOutData(filenameMagneticSusceptibility_X, headerWithBeta, "J1/J2", "magnetic susceptibility in J2", *outDataMagneticSusceptibility_X, N);

        delete outDataMagneticSusceptibility_X;
    }

}