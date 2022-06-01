#include "helpers.h"

/////////////////////////////// bits ///////////////////////////////

void printBits(int s, int N) {
    for (int n = N-1; n >= 0; n--) {
        std::cout << ( (s >> n) & 1 );
    } std::cout << "\n";
}

/////////////////////////////// statistics ///////////////////////////////

std::tuple<double, double> get_mean_and_se(const std::vector<double> &data) {
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / (double) data.size();
    double se = 0.0;
    for (double d : data) {
        se += std::pow(d - mean, 2);
    }
    se /= (double) data.size();
    se = std::sqrt(se);
    return {mean, se};
}

/////////////////////////////// others ///////////////////////////////

// [executable] N J_START J_END J_COUNT CORES SILENT
void validateInput(int &argc, char* argv[], const unsigned int &cpu_cnt, int &N, int &SIZE, double &J_START, double &J_END,
                   int &J_COUNT, double &T_START, double &T_END, int &T_COUNT, bool &silent, int &cores, bool &plotsIn3D,
                   bool skipSilent, const double &J1, const double &J2, bool &noX) {

    if (argc >= 2) {
        std::string DDD = "3D";
        if (argv[1] != DDD) {
            goto no3D;
        }
        std::cout << "3D plots:\n";
        plotsIn3D = true;
        if ( (std::stoi(argv[2]) % 2 == 0) && (std::stoi(argv[2]) >= 6) && (std::stoi(argv[2]) <= 32) ) {
            N = std::stoi(argv[2]);
        } else {
            std::cout << "invalid chain size, must be even and at least 6, defaulting to " << N << "\n";
        }
        SIZE = (int) pow(2, N);
        std::cout << "N: " << N << "; size: " << SIZE << std::endl;
        if (std::stod(argv[3]) > std::stod(argv[4]) || std::stoi(argv[5]) < 1) {
            std::cout << "range invalid, defaulting...\n";
        } else {
            J_START = std::stod(argv[3]); J_END = std::stod(argv[4]); J_COUNT = std::stoi(argv[5]);
        }
        std::cout << "J_START = " << J_START << ", J_END = " << J_END << " and J_COUNT = " << J_COUNT << " from args\n";

        if (std::stod(argv[6]) > std::stod(argv[7]) || std::stoi(argv[8]) < 1) {
            std::cout << "range invalid, defaulting...\n";
        } else {
            T_START = std::stod(argv[6]); T_END = std::stod(argv[7]); T_COUNT = std::stoi(argv[8]);
        }
        std::cout << "T_START = " << T_START << ", T_END = " << T_END << " and T_COUNT = " << T_COUNT << " from args\n";

        if (argc >= 9) {
            int crs = std::stoi(argv[9]);
            if (crs > 0 && crs <= cpu_cnt) {
                cores = crs;
                std::cout << "using " << cores << " cores\n";
            } else {
                std::cout << "defaulting to using all (" << cores << ") cores\n";
            }
        }

        if (argc >= 10) {
            std::string s1 = "noX";
            std::string s2 = argv[10];
//            std::cout << s2 << "\n";
            if (s1 == s2) {
                noX = true;
            }
        }

        if (skipSilent) {
            return;
        }

        if (argc >= 11) {
            std::string s1 = "silent";
            std::string s2 = argv[11];
            std::cout << s2 << "\n";
            if (s1 == s2) {
                silent = true;
            }
        }

        goto silentEntry;
    }

    no3D:

    std::cout << "2D plots:\n";
    if (argc >= 2) {
        if ( (std::stoi(argv[1]) % 2 == 0) && (std::stoi(argv[1]) >= 6) && (std::stoi(argv[1]) <= 32) ) {
            N = std::stoi(argv[1]);
        } else {
            std::cout << "invalid chain size, must be even and at least 6, defaulting to " << N << "\n";
        }
    }
    SIZE = (int) pow(2, N);
    std::cout << "N: " << N << "; size: " << SIZE << std::endl;
    if (argc >= 5) {
        std::cout << "range given: ";
        if (std::stod(argv[2]) > std::stod(argv[3]) || std::stoi(argv[4]) < 1) {
//            std::cout << argv[2] << "\t" << argv[3] << "\t" << argv[4] << "\n";
            std::cout << "range invalid, defaulting...\n";
        } else {
            J_START = std::stod(argv[2]); J_END = std::stod(argv[3]); J_COUNT = std::stoi(argv[4]);
        }
        std::cout << "START = " << J_START << ", END = " << J_END << " and COUNT = " << J_COUNT << " from args\n";
    } else {
        std::cout << "no range given: ";
        std::cout << "START = " << J_START << ", END = " << J_END << " and COUNT = " << J_COUNT << " from default\n";
    }

    std::cout << "default J1 and J2 for beta plot: J1 = " << J1 << " and J2 = " << J2 << " (currently unchangeable)\n";

    if (argc >= 6) {
        int crs = std::stoi(argv[5]);
        if (crs > 0 && crs <= cpu_cnt) {
            cores = crs;
            std::cout << "using " << cores << " cores\n";
        } else {
            std::cout << "defaulting to using all (" << cores << ") cores\n";
        }
    }

    T_START = J_START;
    T_END = J_END;
    T_COUNT = J_COUNT;

    if (argc >= 7) {
        std::string s1 = "noX";
        std::string s2 = argv[6];
        if (s1 == s2) {
            noX = true;
        }
    }

    if (skipSilent) {
        return;
    }

    if (argc >= 8) {
        std::string s1 = "silent";
        std::string s2 = argv[7];
        if (s1 == s2) {
            silent = true;
        }
    }

    silentEntry:

    if (!silent) {
        std::cout << "continue? (y/n):";
        char c;
        std::cin >> c;
        if (c != 'y') {
            exit(1);
//            wrong_N:
//            std::cout << "Enter new N (must be even ans >= 6):";
//            int N_usr;
//            std::cin >> N_usr;
//            if (N_usr >= 6 && N_usr % 2 == 0) {
//                N = N_usr;
//            } else {
//                goto wrong_N;
//            }
//            wrong_JRANGE:
//            std::cout << "Enter new J_START (J1/J2):";
//            double J_START_usr;
//            std::cin >> J_START_usr;
//            std::cout << "Enter new J_END (J1/J2):";
//            double J_END_usr;
//            std::cin >> J_END_usr;
//            std::cout << "Enter new J_COUNT (number of data-points):";
//            int J_COUNT_usr;
//            std::cin >> J_COUNT_usr;
//            if (J_START_usr <= J_END_usr && J_COUNT_usr >= 1) {
//                J_START = J_START_usr;
//                J_END = J_END_usr;
//                J_COUNT = J_COUNT_usr;
//            } else {
//                goto wrong_JRANGE;
//            }
        }
    }
}

std::string formatTime(std::chrono::duration<double> elapsed_seconds) {
    double raw = elapsed_seconds.count();
    int hours = (int) raw / 3600;
    int minutes = (int) (raw - hours * 3600) / 60;
    double seconds = raw - hours * 3600 - minutes *  60;
    std::string ret;
    if (hours > 0) {ret += std::to_string(hours) + " hours " + std::to_string(minutes) + " minutes ";}
    else if (minutes > 0) {ret += std::to_string(minutes) + " minutes ";}
    ret += std::to_string(seconds) + " seconds";
    return ret;
}
