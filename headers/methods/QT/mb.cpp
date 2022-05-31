#include "mb.h"

/////////////////////////////// magnetization blocks ///////////////////////////////

namespace QT::MB {

    std::vector<Eigen::VectorXcd> getVector(const std::vector<matrixType> &matrixBlocks) {

        // random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        std::vector<Eigen::VectorXcd> vectors;

        // fill vector
        for (const matrixType &block : matrixBlocks) {
            int size = (int) block.rows();
            Eigen::VectorXcd v = Eigen::VectorXcd::Zero(size);
            for (int i = 0; i < size; i++) {
                v(i) = std::complex<double>(randNum(gen), randNum(gen));
            }
            v.normalize();
            vectors.push_back(v);
//            std::cout << v << std::endl;
        }

        return vectors;

    }

    std::vector<matrixType> getHamilton(const double &J1, const double &J2, const int &N, const int &SIZE) {

        std::vector<std::vector<int>> states(N + 1);

        int numberOfStates = 0;

        // get states
        for (int s = 0; s < SIZE; s++) {
            int m = ED::bitSum(s, N);
            states.at(m).emplace_back(s);
        }

        std::vector<matrixType> matList;

        for (int m = 0; m <= N; m++) {
            if (states.at(m).empty()) {continue;}
            Eigen::MatrixXd Mtrx = fillHamiltonBlock(J1, J2, states.at(m), N);
            matList.emplace_back(Mtrx.sparseView());
        }

        matList.shrink_to_fit();
        return matList;

    }

    Eigen::MatrixXd fillHamiltonBlock(const double &J1, const double &J2, const std::vector<int> &states, const int &N) {

        const int statesCount = (int) states.size();
        Eigen::MatrixXd hamiltonBlock = Eigen::MatrixXd::Zero(statesCount,statesCount);

        for (int a = 0; a < statesCount; a++) {
            int s = states.at(a);
            for (int n = 0; n < N / 2; n++) {
                // declaring indices
                int j_0 = 2 * n;
                int j_1 = (j_0 + 1) % N;
                int j_2 = (j_0 + 2) % N;
                // applying H to state s
                if (((s >> j_0) & 1) == ((s >> j_2) & 1)) {
                    hamiltonBlock(a,a) += 0.25 * J1;
                } else {
                    hamiltonBlock(a,a) -=0.25 * J1;
                    int d = s ^ (1 << j_0) ^ (1 << j_2);
                    int b = ED::findState(states, d);
                    hamiltonBlock(a,b) = 0.5 * J1;
                }
                if (((s >> j_0) & 1) == ((s >> j_1) & 1)) {
                    hamiltonBlock(a,a) += 0.25 * J2;
                } else {
                    hamiltonBlock(a,a) -= 0.25 * J2;
                    int d = s ^ (1 << j_0) ^ (1 << j_1);
                    int b = ED::findState(states, d);
                    hamiltonBlock(a,b) = 0.5 * J2;
                }
                if (((s >> j_1) & 1) == ((s >> j_2) & 1)) {
                    hamiltonBlock(a,a) += 0.25 * J2;
                } else {
                    hamiltonBlock(a,a) -= 0.25 * J2;
                    int d = s ^ (1 << j_1) ^ (1 << j_2);
                    int b = ED::findState(states, d);
                    hamiltonBlock(a,b) = 0.5 * J2;
                }
            }
        }
        return hamiltonBlock;
    }

    std::vector<matrixType> getS2(const double &J1, const double &J2, const int &N, const int &SIZE) {

        std::vector<std::vector<int>> states(N + 1);

        int numberOfStates = 0;

        // get states
        for (int s = 0; s < SIZE; s++) {
            int m = ED::bitSum(s, N);
            states.at(m).emplace_back(s);
        }

        std::vector<matrixType> matList;

        for (int m = 0; m <= N; m++) {
            if (states.at(m).empty()) {continue;}
            Eigen::MatrixXd Mtrx = fillS2Block(N, states.at(m));
            matList.emplace_back(Mtrx.sparseView());
        }

        matList.shrink_to_fit();
        return matList;

    }

    Eigen::MatrixXd fillS2Block(const int &N, const std::vector<int> &states) {

        int size = (int) states.size();
        Eigen::MatrixXd S2 = 0.75 * (double) N * Eigen::MatrixXd::Identity(size, size);
        for (int k = 0; k < size; k++) {
            for (int j = 0; j < N; j++) {
                int s = states.at(k);
                for (int i = 0; i < j; i++) {
                    if (((s >> i) & 1) == ((s >> j) & 1)) {
                        S2(k, k) += 0.5;
                    } else {
                        S2(k, k) -= 0.5;
                        int d = s ^ (1 << i) ^ (1 << j);
                        int b = ED::findState(states, d);
                        S2(k, b) = 1.0;
                    }
                }
            }
        }
        return S2;
    }

    std::vector<double> rungeKutta4_X(const double &start, const double &end, const double &step, const int &N, const std::vector<matrixType> &H_List, const std::vector<matrixType> &S2_List) {

        std::vector<double> outData;
        std::vector<Eigen::VectorXcd> vec = getVector(H_List);
        int blockCount = (int) H_List.size();

        double norm = 0.0;
        for (int i = 0; i < blockCount; i++) {
            norm += std::pow(vec.at(i).norm(), 2);
        } norm = std::sqrt(norm);
        for (int i = 0; i < blockCount; i++) {
            vec.at(i) = vec.at(i) / norm;
        }

        double beta = start - step;
        while (beta <= end) {

            beta += step;

            norm = 0.0;

#pragma omp parallel for default(none) shared(blockCount, vec, H_List, step, norm)
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = hlp::rungeKutta4Block(vec.at(i), H_List.at(i), step);
                double normnt = std::pow(vec.at(i).norm(), 2);
#pragma omp critical
                norm += normnt;
            }

            norm = std::sqrt(norm);
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = vec.at(i) / norm;
            }

//            norm = 0.0;
//            for (int i = 0; i < blockCount; i++) {
//                norm += std::pow(vec.at(i).norm(), 2);
//            }

            std::vector<double> vec_S2_vec_List;

            double C = 0.0;
#pragma omp parallel for default(none) shared(blockCount, H_List, S2_List, vec, vec_S2_vec_List)
            for (int i = 0; i < blockCount; i++) {
                const matrixType &S2 = S2_List.at(i);
                const Eigen::VectorXcd &v = vec.at(i);
                double v_S2_v = std::real((v.adjoint() * S2 * v)(0,0));
#pragma omp critical
                vec_S2_vec_List.emplace_back(v_S2_v);
            }

            double vec_S2_vec = std::accumulate(vec_S2_vec_List.begin(), vec_S2_vec_List.end(), 0.0);

            C += beta * vec_S2_vec / 3.0 / (double) N;
            outData.emplace_back(C);

        }

        outData.shrink_to_fit();
        return outData;

    }

    void start_calculation_X_J_const(const double &start, const double &end, const double &step,
                                     const double &J1, const double &J2, const int &N, const int &SIZE, const int &SAMPLES) {

        auto start_timer = std::chrono::steady_clock::now();

        std::cout << "\n" << "X(T), J = const, QT, magnetization blocks ..." << std::endl;


        std::vector<matrixType> H_List = getHamilton(J1, J2, N, SIZE);
        std::vector<matrixType> S2_List = getS2(J1, J2, N, SIZE);

        typedef std::vector<std::tuple<double, double>> dataVectorType;
        std::cout << 1 << "/" << SAMPLES << "\n";
        std::vector<double> outData = rungeKutta4_X(start, end, step, N, H_List, S2_List);

        for (int s = 2; s <= SAMPLES; s++) {
            std::cout << s << "/" << SAMPLES << "\n";
            std::vector<double> rawData = rungeKutta4_X(start, end, step, N, H_List, S2_List);
            for (int i = 0; i < rawData.size(); i++) {
                outData.at(i) += rawData.at(i);
            }
        }

        outData.shrink_to_fit();

        for (double &dataPoint : outData) {
            dataPoint = dataPoint / (double) SAMPLES;
        }

        std::vector<double> betaData;
        double beta = start - step;
        while (beta <= end) {
            beta += step;
            betaData.emplace_back(beta);
        } betaData.shrink_to_fit();

//        std::cout << "sizes: " << betaData.size() << "\t" << outData.size() << "\n";

        hlp::saveOutData("data_susceptibility_J_const_QT.txt", "QT, MS für N = " + std::to_string(N)
                                                               + " mit " + std::to_string(SAMPLES) + " Samples",
                         "T in J2 / kb", "C in J2", betaData, outData, N);

        auto end_timer = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
        std::cout << "calculations done; this took: " << formatTime(elapsed_seconds) << "\n";

    }

}
