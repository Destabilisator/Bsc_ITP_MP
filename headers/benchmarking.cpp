#include "benchmarking.h"

namespace bench {

    /////////////////////////////// benchmarking ///////////////////////////////

    void bench_ED_QT_SG(int N_start, int N_end) {
        double J_START = 0.0, J_END = 2.0;
        int J_COUNT = 50;
        double T_START = 0.0, T_END = 50.0;
        for (int N = N_start; N <= N_end; N += 2) {
            int SIZE = (int) std::pow(2, N);
            for (int cores: {1, 2, 5, 10}) {
                for (double stepsize: {0.1, 0.01}) {
                    int T_COUNT = (int) ((T_END - T_START) / stepsize);
                    for (int SAMPLES: {1, 2, 3}) {
                        std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << "), stepsize = " << stepsize
                                  << ", SAMPLES = " << SAMPLES << ", cores = " << cores;
                        // QT
                        auto start_timer_QT = std::chrono::steady_clock::now();
                        QT::MS::start_calc_spin_gap(J_START, J_END, J_COUNT, T_START, T_END, stepsize, N, SIZE, SAMPLES,
                                                    cores);
                        auto end_timer_QT = std::chrono::steady_clock::now();
                        std::chrono::duration<double> elapsed_seconds_QT = end_timer_QT - start_timer_QT;
                        save_bench_val(
                                "QT_SG_step_" + std::to_string(stepsize) + "_SAMPLES_" + std::to_string(SAMPLES) +
                                "_cores_" + std::to_string(cores) + ".txt",
                                std::to_string(N) + "\t" + std::to_string(elapsed_seconds_QT.count()));
                    }
                    // ED with fit
                    std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << "), stepsize = " << stepsize << ", cores = " << cores;
                    auto start_timer_ED_MJ = std::chrono::steady_clock::now();
                    ED::multi::startSusceptibilityMultiJ(J_START, J_END, J_COUNT, T_START, T_END, T_COUNT, N, SIZE,
                                                         cores);
                    auto end_timer_ED_MJ = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed_seconds_ED_MJ = end_timer_ED_MJ - start_timer_ED_MJ;
                    save_bench_val(
                            "ED_MJ_step_" + std::to_string(stepsize) + "_cores_" + std::to_string(cores) + ".txt",
                            std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_MJ.count()));
                }
                // ED from EV
                std::cout << "\nBENCHMARKING: N = " << N << " (" << SIZE << ")" << ", cores = " << cores;
                auto start_timer_ED_SG = std::chrono::steady_clock::now();
                ED::multi::start_SpinGap(J_COUNT, J_START, J_END, cores, N, SIZE);
                auto end_timer_ED_SG = std::chrono::steady_clock::now();
                std::chrono::duration<double> elapsed_seconds_ED_SG = end_timer_ED_SG - start_timer_ED_SG;
                save_bench_val("ED_SG_cores_" + std::to_string(cores) + ".txt",
                               std::to_string(N) + "\t" + std::to_string(elapsed_seconds_ED_SG.count()));
            }
        }
    }

    void save_bench_val(const std::string &filename, const std::string &content) {
        std::ofstream outfile;
        outfile.open("./results/benchmarking/data/" + filename, std::ios_base::app);
//        if (outfile.fail()) {
//            outfile.open("./benchmarking/" + filename);
//        }
        outfile << content << "\n";
        outfile.close();
    }


}
