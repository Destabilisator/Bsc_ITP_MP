#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <thread>
#include <mutex>
#include <ctime>
#include <complex>
//#include <Eigen/Eigenvalues>
//#include <Eigen/SparseCore>

std::mutex coutMutex, nextJMutex;

int CURRENT = 1;

int N = 10;
int size = std::pow(2, N);

int PROGRESSBAR_SEGMENTS = 20;

void printBits(int s, int N) {
    for (int n = N-1; n >= 0; n--) {
        std::cout << ( (s >> n) & 1 );
    } std::cout << "\n";
}

int bitSum(int s, int N) {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        sum += ((s >> i) & 1 );
    } return sum;
}

void fillStatesVector(std::vector<int> *states, int m, int N, int size) {
    for (int s = 0; s <= size - 1; s++) {
        if (bitSum(s, N) == m) {
            states->push_back(s);
        }
    } states->shrink_to_fit();
}

void fillStatesList(std::list<int> *states, int m, int N, int size) {
    for (int s = 0; s <= size - 1; s++) {
        if (bitSum(s, N) == m) {
            states->push_back(s);
        }
    }
}

void get_DeltaE_CT_const(double J, int pos, const int &COUNT, const double &START, const double &END, const int &N, const int &SIZE) {

    // progressbar init
    coutMutex.lock();
    std::cout << "\r[";
    for (int _ = 0; _ < PROGRESSBAR_SEGMENTS; _++) {
        std::cout << ".";
    } std::cout << "] " << int(0.0) << "% J1/J2 = " << START << " (" << 0 << "/" << COUNT << ")     ";
    std::cout.flush();
    coutMutex.unlock();

    while (true) {

        // progressbar
        nextJMutex.lock();
        int prg = std::min(CURRENT, COUNT);
        int p = (int) ( (float) prg / (float) COUNT * (float) PROGRESSBAR_SEGMENTS);
        //coutMutex.lock();
        std::cout << "\r[";
        for (int _ = 0; _ < p; _++) {
            std::cout << "#";
        } for (int _ = p; _ < PROGRESSBAR_SEGMENTS; _++) {
            std::cout << ".";
        } std::cout << "] " << int( (float) prg / (float) COUNT * 100.0 ) << "% J1/J2 = " << J << " (" << prg << "/" << COUNT << ")     ";
        std::cout.flush();
        //coutMutex.unlock();

        // write data
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


void start_DeltaE_CT_const(const int &COUNT, const double &START, const double &END, const unsigned int &cpu_cnt,
                           int &cores, const int &N, const int &SIZE) {

    auto start = std::chrono::steady_clock::now();

    std::cout << "\n" << "Delta E and C (T=const): calculating:...";

    if (COUNT < cores) {
        cores = COUNT;
    }
    std::cout << " (" << cores << ") cores\n";

    std::thread Threads[cores];

    CURRENT = 1 + cores;

    for (int i = 0; i < cores; i++) {
        Threads[i] = std::thread(get_DeltaE_CT_const, START + (END - START) * i / COUNT, i + 1, COUNT, START, END, N, SIZE);
    }

    for (int i = 0; i < cores; i++) {
        Threads[i].join();
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "\n" << "calculations done; this took: " << elapsed_seconds.count() << " seconds\n\n";

    std::cout << "\n";

}

void fillStateBench() {
    std::ofstream file;
    try {
        file.open("benchmarkResults.txt");
        file << "N\t" << "fillstatesVector | " << "fillstatesList | " << "assign states vector | "
             << "assign states vector list" << "\n\n";

        for (int N = 0; N <= 32; N++) {

            std::cout << "benchmarking N = " << N << std::endl;
            size = pow(2, N);
            file << N << "\t";

            // fill states
            const clock_t begin_time_BENCH_1 = clock();
            for (int m = 0; m <= N; m++) {
                auto *states = new std::vector<int>;
                fillStatesVector(states, m, N, size);
            }

            double time_BENCH_1 = float(clock() - begin_time_BENCH_1) / CLOCKS_PER_SEC;
            std::cout << "benchmark 1 (filling states vector) done; this took: " << time_BENCH_1 << " seconds\n";
            file << time_BENCH_1 << "\t";


            // fill states
            const clock_t begin_time_BENCH_1_1 = clock();
            for (int m = 0; m <= N; m++) {
                auto *states = new std::list<int>;
                fillStatesList(states, m, N, size);
            }

            double time_BENCH_1_1 = float(clock() - begin_time_BENCH_1_1) / CLOCKS_PER_SEC;
            std::cout << "benchmark 1_1 (filling states list) done; this took: " << time_BENCH_1_1 << " seconds\n";
            file << time_BENCH_1_1 << "\t";


            // assigning states
            const clock_t begin_time_BENCH_2 = clock();
            auto *states = new std::vector<std::vector<int>>(N + 1);
            for (int s = 0; s < size; s++) {
                states->at(bitSum(s, N)).push_back(s);
            }

            double time_BENCH_2 = float(clock() - begin_time_BENCH_2) / CLOCKS_PER_SEC;
            std::cout << "benchmark 2 (assigning states vector) done; this took: " << time_BENCH_2 << " seconds\n";
            file << time_BENCH_2 << "\t";


            // assigning states
            const clock_t begin_time_BENCH_2_1 = clock();
            auto *statesList = new std::vector<std::list<int>>(N + 1);
            for (int s = 0; s < size; s++) {
                statesList->at(bitSum(s, N)).push_back(s);
            }

            double time_BENCH_2_1 = float(clock() - begin_time_BENCH_2_1) / CLOCKS_PER_SEC;
            std::cout << "benchmark 2_1 (assigning states vector list) done; this took: " << time_BENCH_2_1
                      << " seconds\n";
            file << time_BENCH_2_1 << "\n";

            std::cout << std::endl;

        }
    } catch (...) {
        std::cout << "interrupted, closing file...\n";
        file.close();
    }

}



int main() {

    int COUNT = 100000;
    double START = 0.0;
    double END = 2.0;
    unsigned int cpu_cnt = 16;
    int cores = 16;


    start_DeltaE_CT_const(COUNT, START, END, cpu_cnt, cores, N, size);

    return 0;
}
