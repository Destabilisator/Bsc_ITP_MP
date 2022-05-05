#include "main.h"

/////////////////////////////// MULTITHREADING ///////////////////////////////

#ifdef multiCalc
void thread_function(double J, int J_pos, std::vector<std::tuple<double, double>> *outDataDeltaE,
                     double beta, std::vector<std::tuple<double, double>> *outDataSpecificHeat_C,
                     std::vector<std::tuple<double, double>> *outDataMagneticSusceptibility_X) {

    while (true) {

        int p = (int) ( (float) J_pos / (float) J_COUNT * (float) PROGRESSBAR_SEGMENTS);
        coutMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < p; _++) {
            std::cout << "#";
        } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) J_pos / (float) J_COUNT * 100.0 ) << "% J1/J2 = " << J << " (" << J_pos << "/" << J_COUNT << ")     ";
        coutMutex.unlock();

        auto *eiVals = new std::vector<std::complex<double>>;
        auto *matrixBlocks = new std::vector<Eigen::MatrixXcd>;
        //auto *matrixBlocks = new std::vector<Eigen::MatrixXd>;

        //naiv::getEiVals(J, 1.0, eiVals, N, SIZE);
        //magnetizationBlocks::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);
        momentumStates::getEiVals(J, 1.0, eiVals, matrixBlocks, N, SIZE);

        // sort eigenvalues
        eiVals->shrink_to_fit();
        std::sort(eiVals->begin(), eiVals->end(), [](const std::complex<double> &c1, const std::complex<double> &c2) {
            return std::real(c1) < std::real(c2);
        });

        // Delta E
        double E0 = std::real(eiVals->at(0));
        double E1 = std::real(eiVals->at(1));

        ///// magnetic susceptibility ////////// magnetic susceptibility ////////// magnetic susceptibility /////
        double magneticSusceptibility_X = 0.0;
        outDataMagneticSusceptibility_X->push_back({J, magneticSusceptibility_X});

        // write data
        nextJMutex.lock();
            outDataDeltaE->push_back({J, E1 - E0});
            outDataSpecificHeat_C->push_back({J, getSpecificHeat(beta, *eiVals, N)});
            //outDataMagneticSusceptibility_X->push_back({J, magneticSusceptibility_X});
            J_pos = J_CURRENT;
            J_CURRENT++;
        nextJMutex.unlock();

        if (J_pos > J_COUNT) {
            delete eiVals;
            delete matrixBlocks;
            break;
        } else {
            J = J_START + (J_END-J_START)*J_pos/J_COUNT;
        }
        // clean up
        eiVals->clear();
        matrixBlocks->clear();
    }

}
#endif

/////////////////////////////// MAIN ///////////////////////////////

int main(int argc, char* argv[]) {

    bool silent = false;
    int cores = (int) cpu_cnt;

    validateInput(argc, argv, &N, &SIZE, &J_START, &J_END, &J_COUNT, &cpu_cnt, &silent, &cores, &J1, &J2);

    // syncing BETA and J ranges
    BETA_START = J_START;
    BETA_END = J_END;
    BETA_COUNT = J_COUNT;

/////////////////////////////// MULTITHREADING ///////////////////////////////

#ifdef multiCalc
    const clock_t begin_time = clock();

    std::cout << "\ncalculating:..." << std::endl;

    //auto *outDataDeltaE = new std::vector<std::tuple<double, std::complex<double>>>;
    auto *outDataDeltaE = new std::vector<std::tuple<double, double>>;
    auto *outDataSpecificHeat_C = new std::vector<std::tuple<double, double>>;
    auto *outDataMagneticSusceptibility_X = new std::vector<std::tuple<double, double>>;

    if (J_COUNT < cpu_cnt) {
        cores = J_COUNT;
    }
    std::thread Threads[cores];

    J_CURRENT += cores;

    for (int i = 0; i < cores; i++) {
        Threads[i] = std::thread(thread_function, J_START + (J_END - J_START) * i / J_COUNT, i + 1, outDataDeltaE, BETA,
                                 outDataSpecificHeat_C, outDataMagneticSusceptibility_X);
    }

    for (int i = 0; i < cores; i++) {
        Threads[i].join();
    }

    auto time = float(clock () - begin_time) /  CLOCKS_PER_SEC;
    std::cout << "\n" << "calculations done; this took: " << time << " seconds\n\n";

    // sort data-points
    std::sort(outDataDeltaE->begin(), outDataDeltaE->end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
        return std::get<0>(a) < std::get<0>(b);
    });
    std::sort(outDataSpecificHeat_C->begin(), outDataSpecificHeat_C->end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
        return std::get<0>(a) < std::get<0>(b);
    });
//    std::sort(outDataMagneticSusceptibility_X->begin(), outDataMagneticSusceptibility_X->end(), [](const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
//        return std::get<0>(a) < std::get<0>(b);
//    });

    std::string filenameDeltaE = "data_delta_E.txt";
    std::string filenameSpecificHeat_C = "data_specific_heat.txt";
    //std::string filenameMagneticSusceptibility_X = "data_magnetic_susceptibility.txt";
    std::string header = "N: " + std::to_string(N) + "\n"
                         + "J1/J2 START: " + std::to_string(J_START) + "\n"
                         + "J1/J2 END: " + std::to_string(J_END) + "\n"
                         + "data-points: " + std::to_string(J_COUNT) + "\n"
                         + "calculation time with " + std::to_string(cores) + " threads: " + std::to_string(time) + " seconds";

    std::string headerWithBeta = "BETA = " + std::to_string(BETA) + "\n" + header;

    saveOutData(filenameDeltaE, header, "J1/J2", "Delta E in J2", *outDataDeltaE, N);
    saveOutData(filenameSpecificHeat_C, headerWithBeta, "J1/J2", "specific heat in J2", *outDataSpecificHeat_C, N);
    //saveOutData(filenameMagneticSusceptibility_X, headerWithBeta, "J1/J2", "magnetic susceptibility in J2", *outDataMagneticSusceptibility_X);

#endif

/////////////////////////////// calculate quantities ///////////////////////////////

    momentumStates::startSpecificHeat(J1, J2, N, SIZE, BETA_START, BETA_END, BETA_COUNT, cores);
    naiv::start(J1, J2, N, SIZE, BETA_START, BETA_END, BETA_COUNT, cores);
    magnetizationBlocks::startSusceptibility(J1, J2, N, SIZE, BETA_START, BETA_END, BETA_COUNT, cores);

/////////////////////////////// naiver Ansatz ///////////////////////////////

#ifdef naiverAnsatz
    naiv::start(J1, J2, N, SIZE, BETA_START, BETA_END, BETA_COUNT, cores);
#endif

/////////////////////////////// fixed magnetization blocks ///////////////////////////////

#ifdef magnetizationBlocksAnsatz
    magnetizationBlocks::start(J1, J2, N, SIZE);
#endif

/////////////////////////////// momentum states ///////////////////////////////

#ifdef momentumStateAnsatz
    momentumStates::start(J1, J2, N, SIZE);
#endif

/////////////////////////////// parity states (unfinished) ///////////////////////////////

#ifdef parityStateAnsatz
    parityStates::start(J1, J2, N, SIZE);
#endif

    return 0;
}
