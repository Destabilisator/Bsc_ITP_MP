#include <iostream>
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

void fillStates(std::vector<int> *states, int m, int N, int size) {
    for (int s = 0; s <= size - 1; s++) {
        if (bitSum(s, N) == m) {
            states->push_back(s);
        }
    } states->shrink_to_fit();
}

int main() {

    for (int N = 0; N <= 32 ; N ++) {
        std::cout << "benchmarking N = " << N << std::endl;
        size = pow(2, N);
        // fill states
        const clock_t begin_time_BENCH_1 = clock();
        for (int m = 0; m <= N; m++) {
            auto *states = new std::vector<int>;
            fillStates(states, m, N, size);
        }

        auto time_BENCH_1 = float( clock () - begin_time_BENCH_1 ) /  CLOCKS_PER_SEC;
        std::cout << "benchmark 1 (filling states) done; this took: " << time_BENCH_1 << " seconds\n";


        // assigning states
        const clock_t begin_time_BENCH_2 = clock();
        auto *states = new std::vector<std::vector<int>>(N+1);
        for (int s = 0; s < size; s++) {
            states->at(bitSum(s, N)).push_back(s);
        }

        auto time_BENCH_2 = float( clock () - begin_time_BENCH_2 ) /  CLOCKS_PER_SEC;
        std::cout << "benchmark 2 (assigning states) done; this took: " << time_BENCH_2 << " seconds\n";

        std::cout << std::endl;

        /*
        for (std::vector<int> vec : *states) {
            for (int s : vec) {
                printBits(s, N);
            }
        }*/
    }

    return 0;
}
