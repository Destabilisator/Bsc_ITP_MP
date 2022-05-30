#include "multi.h"
/*
/////////////////////////////// multi-threading ///////////////////////////////

namespace QT::multi {

    /////////////////////////////// C(T) ///////////////////////////////

    void start_C_J_const(const double &START, const double &END, const double &STEP, int &cores,
                         const double &J, const int &N, const int &SIZE) {

        auto start = std::chrono::steady_clock::now();

        std::cout << "\n" << "Delta E and C (T=const): calculating:...";

        std::vector<indexStateVectorType> indexList = MS::getIndexAndStates(N, SIZE);
        std::vector<std::vector<std::tuple<double, double>>> rawData;

        if ((int) indexList.size() < cores) {
            cores = (int) indexList.size();
        }
        std::cout << " (" << cores << ") cores";

        std::thread Threads[cores];

        NEXT = 0 + cores;

        std::cout << ", momentum states\n";
        for (int i = 0; i < cores; i++) {
            Threads[i] = std::thread(get_C_J_const, J, i , &rawData, indexList, START, END, STEP, N, SIZE);
        }

        for (int i = 0; i < cores; i++) {
            Threads[i].join();
        }

        std::vector<std::tuple<double, double>> outData;

        for (int i = 0; i < rawData.at(0).size(); i++) {
            double temp = std::get<0>(rawData.at(0).at(i));
            double C = 0.0;
            for (const std::vector<std::tuple<double, double>> &dat : rawData) {
                C += std::get<1>(dat.at(i));
            }
            C = C / (double) N;
//            std::cout << temp << " " << C << "\n";
            outData.emplace_back(temp, C);
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "\n" << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n\n";

        // sort data-points
        std::sort(outData.begin(), outData.end(), [](
                const std::tuple<double, double> &a, const std::tuple<double, double> &b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::string filename = "data_specific_heat_J_const_QT_multi.txt";
        std::string header = "N: " + std::to_string(N) + "\n"
                             + "T START: " + std::to_string(START) + "\n"
                             + "T END: " + std::to_string(END) + "\n"
                             + "step size: " + std::to_string(STEP) + "\n"
                             + "calculation time with " + std::to_string(cores) + " threads: " +
                             std::to_string(elapsed_seconds.count()) + " seconds";

        std::string headerWithJ = "J = " + std::to_string(J) + "\n" + header;

        hlp::saveOutData(filename, headerWithJ, "J1/J2", "specific heat in J2", outData, N);

    }

    void get_C_J_const(double J, int pos, std::vector<std::vector<std::tuple<double, double>>> *outData,
                       const std::vector<indexStateVectorType> &indexList, const double &START,
                       const double &END, const double &STEP, const int &N, const int &SIZE) {

        const int COUNT = (int) indexList.size() - 1;

        // progressbar init
        nextQTMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << 0 <<  "% block = " << 0 << "/" << COUNT << "     ";
        std::cout.flush();
        nextQTMutex.unlock();

        while (pos < (int) indexList.size()) {

            std::vector<std::tuple<double, double>> dat = MS::rungeKutta4_C(START, END, STEP, J, 1.0, N, indexList.at(pos));

            // progressbar
            nextQTMutex.lock();
            int prg = std::min({NEXT, COUNT});;
            int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% block = " << prg << "/" << COUNT << "     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            outData->push_back(dat);
            pos = NEXT;
            NEXT++;
            nextQTMutex.unlock();

        }

    }

}*/