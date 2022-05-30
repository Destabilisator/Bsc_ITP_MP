#include "multi.h"

/////////////////////////////// multi-threading ///////////////////////////////
/*
namespace QT::multi {

    /////////////////////////////// C(T) ///////////////////////////////

    void start_C_3D(const double &J_START, const double &J_END, const int &J_COUNT,
                    const double &T_START, const double &T_END, const double &T_STEP,
                    int &cores, const int &N, const int &SIZE) {

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
            double J = (double) J_START + (double) (J_END - J_START) * (double) i / (double) J_COUNT;
            Threads[i] = std::thread(get_C_J_const, J_START, J_END, J_COUNT, T_START, T_END, T_STEP, J, i, N, SIZE);
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

    }

    void get_C_J_3D(const double &J_START, const double &J_END, const int &J_COUNT,
                    const double &T_START, const double &T_END, const double &T_STEP,
                    double J, int pos, const int &N, const int &SIZE) {

        // progressbar init
        nextQTMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int(0.0) << "% J1/J2 = " << J_START << " (" << 0 << "/" << J_COUNT << ")     ";
        std::cout.flush();
        nextQTMutex.unlock();

        while (pos < J_COUNT) {

            std::vector<matrixType> matrixList = QT::MS::getHamilton(J, 1.0, N, SIZE);

            std::vector<std::tuple<double, double>> data = QT::MS::rungeKutta4_C(T_START, T_END, T_STEP, J, 1.0, N, SIZE, matrixList);

            // progressbar
            nextQTMutex.lock();
            int prg = std::min({NEXT, J_COUNT});;
            int p = (int) ( (float) prg / (float) J_COUNT * (float) PROGRESSBAR_SEGMENTS);
            //coutMutex.lock();
            std::cout << "\r[";
            for (int _ = 0; _ < p; _++) {
                std::cout << "#";
            } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
                std::cout << ".";
            } std::cout << "] " << int( (float) prg / (float) J_COUNT * 100.0 ) << "% block = " << prg << "/" << J_COUNT << "     ";
            std::cout.flush();
            //coutMutex.unlock();

            // write data
            pos = NEXT;
            NEXT++;
            nextQTMutex.unlock();

            J = (double) J_START + (double) (J_END - J_START) * (double) pos / (double) J_COUNT;

        }

    }

}*/