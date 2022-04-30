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

int N, size;

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

int main() {

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

    return 0;
}
