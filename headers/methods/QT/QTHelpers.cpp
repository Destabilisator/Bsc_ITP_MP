#include "QTHelpers.h"

namespace QT::hlp {

    ///// random vectors /////

    Eigen::VectorXcd getVector(int size) {

        // random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        // fill vector
        Eigen::VectorXcd v = Eigen::VectorXcd::Zero(size);
        for (int i = 0; i < size; i++) {
            v(i) = std::complex<double>(randNum(gen), randNum(gen));
        }
        v.normalize();
        return v;

    }

/*
    std::vector<Eigen::VectorXcd> getVector(const int &N, const int &SIZE, const std::vector<matrixDataMomentumType> &matrixBlocks) {

        // random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        std::vector<Eigen::VectorXcd> vectors;

        // fill vector
        for (const matrixDataMomentumType &data : matrixBlocks) {
            int size = (int) std::get<2>(data).size();
            Eigen::VectorXcd v = Eigen::VectorXcd::Zero(size);
            for (int i = 0; i < size; i++) {
                v(i) = std::complex<double>(randNum(gen), randNum(gen));
            }
            v.normalize();
            vectors.push_back(v);
            std::cout << v << std::endl;
        }

        return vectors;

    }
*/

    std::vector<Eigen::VectorXcd> getVector(const std::vector<matrixTypeComplex> &matrixBlocks) {

        // random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        std::vector<Eigen::VectorXcd> vectors;

        // fill vector
        for (const matrixTypeComplex &block : matrixBlocks) {
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

    std::vector<Eigen::VectorXd> getVector(const std::vector<matrixTypeReal> &matrixBlocks) {

        // random number generation
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> randNum{0,1};

        std::vector<Eigen::VectorXd> vectors;

        // fill vector
        for (const matrixTypeReal &block : matrixBlocks) {
            int size = (int) block.rows();
            Eigen::VectorXd v = Eigen::VectorXd::Zero(size);
            for (int i = 0; i < size; i++) {
                v(i) = randNum(gen);
            }
            v.normalize();
            vectors.push_back(v);
//            std::cout << v << std::endl;
        }

        return vectors;

    }

    /////////////////////////////// Runge-Kutta 4 ///////////////////////////////

    Eigen::VectorXcd rungeKutta4Block(const Eigen::VectorXcd &vec, const matrixTypeComplex &H, const double &step) {
        Eigen::VectorXcd k1 = - 0.5 * H * vec;
        Eigen::VectorXcd k2 = - 0.5 * H * (vec + step * 0.5 * k1);
        Eigen::VectorXcd k3 = - 0.5 * H * (vec + step * 0.5 * k2);
        Eigen::VectorXcd k4 = - 0.5 * H * (vec + step * k3);

        return vec + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

    Eigen::VectorXd rungeKutta4Block(const Eigen::VectorXd &vec, const matrixTypeReal &H, const double &step) {
        Eigen::VectorXd k1 = - 0.5 * H * vec;
        Eigen::VectorXd k2 = - 0.5 * H * (vec + step * 0.5 * k1);
        Eigen::VectorXd k3 = - 0.5 * H * (vec + step * 0.5 * k2);
        Eigen::VectorXd k4 = - 0.5 * H * (vec + step * k3);

        return vec + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

    std::vector<Eigen::VectorXcd> normalizedVectorList(const std::vector<Eigen::VectorXcd> &vectors, const int &SIZE) {

        Eigen::VectorXcd v = Eigen::VectorXcd::Zero(SIZE);

        // put together v_ret
        int pos = 0;
        for (const Eigen::VectorXcd &vector : vectors) {
            for (int i = 0; i < vector.rows(); i++) {
                v(pos) = vector(i);
                pos++;
            }
        } v.normalize();

        // put together vectors list
        std::vector<Eigen::VectorXcd> vectorList;
        pos = 0;
        for (const Eigen::VectorXcd &vector : vectors) {
            int sz = (int) vector.size();
            Eigen::VectorXcd v_ls = Eigen::VectorXcd::Zero(sz);
            for (int j = 0; j < sz; j++) {
                v_ls(j) = v(pos);
                pos++;
            }
            vectorList.emplace_back(v);
        }

        return vectorList;

    }

    ///// C /////

    std::vector<double> rungeKutta4_C(const double &start, const double &end, const double &step, const int &N,
                                      const std::vector<matrixTypeComplex> &matrixList) {

        std::vector<double> outData;
        std:std::vector<Eigen::VectorXcd> vec = getVector(matrixList);
        int blockCount = (int) matrixList.size();

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
//#if INNER_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, vec, matrixList, step, norm)
//#endif
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = hlp::rungeKutta4Block(vec.at(i), matrixList.at(i), step);
                double normnt = std::pow(vec.at(i).norm(), 2);
//#pragma omp critical
                norm += normnt;
            }

            norm = std::sqrt(norm);
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = vec.at(i) / norm;
            }

            std::vector<double> vec_H_vec_List;
            std::vector<double> vec_H2_vec_List;

            //double C = 0.0;
//#if INNER_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, matrixList, vec, vec_H_vec_List, vec_H2_vec_List)
//#endif
            for (int i = 0; i < blockCount; i++) {
                const matrixTypeComplex &H = matrixList.at(i);
                const Eigen::VectorXcd &v = vec.at(i);
                double vHv = std::real((v.adjoint() * H * v)(0,0));
                double vH2v = std::real((v.adjoint() * H * H * v)(0,0));
//#pragma omp critical
                vec_H_vec_List.emplace_back(vHv);
//#pragma omp critical
                vec_H2_vec_List.emplace_back(vH2v);
            }

            double vec_H_vec = std::accumulate(vec_H_vec_List.begin(), vec_H_vec_List.end(), 0.0);
            double vec_H2_vec = std::accumulate(vec_H2_vec_List.begin(), vec_H2_vec_List.end(), 0.0);

            double H_diff = std::real(vec_H2_vec - std::pow(vec_H_vec, 2));
            double C = beta * beta * H_diff / (double) N;
//#pragma omp critical
            outData.emplace_back(C);

        }

        outData.shrink_to_fit();
        return outData;

    }

    ///// X /////

    std::vector<double> rungeKutta4_X(const double &start, const double &end, const double &step, const int &N,
                                      const std::vector<matrixTypeComplex> &H_List, const std::vector<matrixTypeComplex> &S2_List) {

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

            // RK4 with vec on H to get new state
//#if INNER_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, vec, H_List, step, norm)
//#endif
//#pragma omp parallel for default(none) shared(blockCount, vec, H_List, step, norm)
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = hlp::rungeKutta4Block(vec.at(i), H_List.at(i), step);
//                double normnt = std::pow(vec.at(i).norm(), 2);
//#pragma omp critical
                norm += std::pow(vec.at(i).norm(), 2);
            }

            // ensure norm(vec) = 1
            norm = std::sqrt(norm);
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = vec.at(i) / norm;
            }

            std::vector<double> vec_S2_vec_List;

//#if INNER_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, S2_List, vec, vec_S2_vec_List)
//#endif
//#pragma omp parallel for default(none) shared(blockCount, S2_List, vec, vec_S2_vec_List)
            for (int i = 0; i < blockCount; i++) {
                const matrixTypeComplex &S2 = S2_List.at(i);
                const Eigen::VectorXcd &v = vec.at(i);
                double vS2v = std::real((v.adjoint() * S2 * v)(0,0));
//#pragma omp critical
                vec_S2_vec_List.emplace_back(vS2v);
            }

            double vec_S2_vec = std::accumulate(vec_S2_vec_List.begin(), vec_S2_vec_List.end(), 0.0);

            double X = beta * vec_S2_vec / 3.0 / (double) N;

            outData.emplace_back(X);

        }

        outData.shrink_to_fit();
        return outData;

    }

    std::vector<Trp> spinMatrixMomentum(const int &N, const int &k, const std::vector<int> &states, const std::vector<int> &R_vals) {

        const int size = (int) states.size();
        std::vector<Trp> S2;
        for (int a = 0; a < size; a++) {
            S2.emplace_back(Trp(a, a, std::complex<double>(0.75 * (double) N, 0.0)));
        }
        for (int a = 0; a < size; a++) {
            int s = states.at(a);
            for (int j = 0; j < N; j++) {
                for (int i = 0; i < j; i++) {
                    if (((s >> i) & 1) == ((s >> j) & 1)) {
                        S2.emplace_back(Trp(a, a, std::complex<double>(0.5, 0.0)));
                    } else {
                        S2.emplace_back(Trp(a, a, std::complex<double>(-0.5, 0.0)));
                        int d = s ^ (1 << i) ^ (1 << j);
                        int r = 0, l = 0;
                        ED::representative(d, &r, &l, N);
                        int b = ED::findState(states, r);
                        if (b >= 0) {
                            std::complex<double> numC(0.0, 4 * PI * (double) k * (double) l / (double) N);
                            S2.emplace_back(Trp(a,b, (std::complex<double>) 1.0 * sqrt((double) R_vals.at(a)
                                                                         / (double) R_vals.at(b)) * std::exp(numC) ));
                        }
                    }
                }
            }
        }
        //std::cout << S2 << std::endl;
        return S2;
    }

    ///// C and X /////

    std::vector<std::tuple<double, double>> rungeKutta4_CX(const double &start, const double &end, const double &step, const int &N,
                                                           const std::vector<matrixTypeComplex> &H_List, const std::vector<matrixTypeComplex> &S2_List) {

        std::vector<std::tuple<double, double>> outData;
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

            // RK4 with vec on H to get new state
//#if INNER_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, vec, H_List, step, norm)
//#endif
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = hlp::rungeKutta4Block(vec.at(i), H_List.at(i), step);
                double normnt = std::pow(vec.at(i).norm(), 2);
//#pragma omp critical
                norm += normnt;
            }

            // ensure norm(vec) = 1
            norm = std::sqrt(norm);
            for (int i = 0; i < blockCount; i++) {
                vec.at(i) = vec.at(i) / norm;
            }

            std::vector<double> vec_H_vec_List;
            std::vector<double> vec_H2_vec_List;
            std::vector<double> vec_S2_vec_List;

//#if INNER_NESTED_THREADS > 1
//#pragma omp parallel for num_threads(INNER_NESTED_THREADS) default(none) shared(blockCount, H_List, S2_List, vec, vec_S2_vec_List, vec_H_vec_List, vec_H2_vec_List)
//#endif
            for (int i = 0; i < blockCount; i++) {
                const matrixTypeComplex &H = H_List.at(i);
                const matrixTypeComplex &S2 = S2_List.at(i);
                const Eigen::VectorXcd &v = vec.at(i);
                double vS2v = std::real((v.adjoint() * S2 * v)(0,0));
                double vHv = std::real((v.adjoint() * H * v)(0,0));
                double vH2v = std::real((v.adjoint() * H * H * v)(0,0));
//#pragma omp critical
                vec_S2_vec_List.emplace_back(vS2v);
//#pragma omp critical
                vec_H_vec_List.emplace_back(vHv);
//#pragma omp critical
                vec_H2_vec_List.emplace_back(vH2v);
            }

            double vec_H_vec = std::accumulate(vec_H_vec_List.begin(), vec_H_vec_List.end(), 0.0);
            double vec_H2_vec = std::accumulate(vec_H2_vec_List.begin(), vec_H2_vec_List.end(), 0.0);
            double H_diff = std::real(vec_H2_vec - std::pow(vec_H_vec, 2));
            double C = beta * beta * H_diff / (double) N;

            double vec_S2_vec = std::accumulate(vec_S2_vec_List.begin(), vec_S2_vec_List.end(), 0.0);
            double X = beta * vec_S2_vec / 3.0 / (double) N;

            outData.emplace_back(C,X);

        }

        outData.shrink_to_fit();
        return outData;

    }

    /////////////////////////////// saving data ///////////////////////////////

    void saveOutData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                     const std::vector<std::tuple<double, double>> &outData, const int &N) {

        std::cout << "saving to file '" << std::to_string(N) << "/data/" << filename << "'... ";// << std::endl;
        std::ofstream file;
        try {
            file.open("./results/" + std::to_string(N) + "/data/" + filename);
            file << header <<"\n\n";
            file << x_lbl << "\t" << y_lbl << "\n";
            for (std::tuple<double, double> data : outData) {
                file << std::get<0>(data) << "\t" << std::get<1>(data) << "\n";
            }
        } catch (...) {
            file.close();
            std::cout << "failed to save to file\n";
        }
        std::cout << "done\n";
        file.close();
    }

    void saveOutData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                     const std::vector<double> &xData, const std::vector<double> &yData, const int &N) {

        std::cout << "saving to file '" << std::to_string(N) << "/data/" << filename << "'... ";// << std::endl;
        std::ofstream file;
        try {
            file.open("./results/" + std::to_string(N) + "/data/" + filename);
            file << header <<"\n\n";
            file << x_lbl << "\t" << y_lbl << "\n";
            for (int i = 0; i < xData.size(); i++) {
                file << xData.at(i) << "\t" << yData.at(i) << "\n";
            }
        } catch (...) {
            file.close();
            std::cout << "failed to save to file\n";
        }
        std::cout << "done\n";
        file.close();
    }

    void saveOutDataSilent(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                           const std::vector<double> &xData, const std::vector<double> &yData, const int &N) {

        std::ofstream file;
        try {
            file.open("./results/" + std::to_string(N) + "/data/" + filename);
            file << header <<"\n\n";
            file << x_lbl << "\t" << y_lbl << "\n";
            for (int i = 0; i < xData.size(); i++) {
                file << xData.at(i) << "\t" << yData.at(i) << "\n";
            }
        } catch (...) {
            file.close();
            std::cout << "failed to save to file\n";
        }
        file.close();
    }

    void saveOutData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                     const std::vector<double> &xData, const std::vector<double> &yData, const std::vector<double> &yErrData, const int &N) {

        std::cout << "saving to file '" << std::to_string(N) << "/data/" << filename << "'... ";// << std::endl;
        std::ofstream file;
        try {
            file.open("./results/" + std::to_string(N) + "/data/" + filename);
            file << header <<"\n\n";
            file << x_lbl << "\t" << y_lbl << "\t" << "yErr" << "\n";
            for (int i = 0; i < xData.size(); i++) {
                file << xData.at(i) << "\t" << yData.at(i) << "\t" << yErrData.at(i) << "\n";
            }
        } catch (...) {
            file.close();
            std::cout << "failed to save to file\n";
        }
        std::cout << "done\n";
        file.close();
    }

    void saveStatisticsData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                            const std::vector<double> &xData, const std::vector<std::vector<double>> &RAWyData,
                            const int &SAMPLES, const double &step, const int &N) {
#ifdef SAVE_WITH_SETS_OF_n_SAMPLES
        // save for different amounts n of random states
        for (int n = 1; n <= SAMPLES / 2; n++) {
            std::vector<std::vector<double>> Y_Data_to_avg;
            // combine data for every i to i+n samples into one avg
            for (int i = 0; i+n < SAMPLES; i += n) {
                std::vector<std::vector<double>> Y_Data_raw;
                // get data to avg
                for (int j = i; j < i + n; j++) {
//                    std:: cout << "n: " << n << ", i: " << i << ", j: " << j << std::endl;
                    Y_Data_raw.emplace_back(RAWyData.at(j));
                }
                std::vector<double> Y_vector;
                // avg
                for(int k = 0; k < Y_Data_raw.at(0).size(); k++) {
                    std::vector<double> Y_vector_raw;
                    for (std::vector<double> Yd : Y_Data_raw) {
                        Y_vector_raw.emplace_back(Yd.at(k));
                    }
                    std::tuple<double, double> mean_se = get_mean_and_se(Y_vector_raw);
                    Y_vector.emplace_back(std::get<0>(mean_se));
                }
                Y_Data_to_avg.emplace_back(Y_vector);
            }
            std::vector<double> Y_Data_save;
            std::vector<double> YErr_Data_save;
            // avg and stdv of avg C with n samples each
            for(int k = 0; k < Y_Data_to_avg.at(0).size(); k++) {
                std::vector<double> Y_vector_raw;
                for (std::vector<double> Yd : Y_Data_to_avg) {
                    Y_vector_raw.emplace_back(Yd.at(k));
                } std::tuple<double, double> mean_se_Y = get_mean_and_se(Y_vector_raw);
                Y_Data_save.emplace_back(std::get<0>(mean_se_Y));
                YErr_Data_save.emplace_back(std::get<1>(mean_se_Y));
            }
            // save data N_n_data_specific_heat_J_const_QT_txt
            std::string flnm = std::to_string(n) + "_" + filename;
            hlp::saveOutData(flnm + ".txt", header, x_lbl, y_lbl, xData, Y_Data_save, YErr_Data_save, N);
#ifdef SAVE_WITH_STEP_SIZE
            hlp::saveOutData("step_size_data/" + flnm + "_step" + std::to_string(step) + ".txt",
                             header, x_lbl, y_lbl, xData, Y_Data_save, YErr_Data_save, N);
#endif
        }

#endif

        #ifdef SAVE_WITH_DATA_FROM_ALL_SAMPLES
        std::vector<double> Y_Data;
        std::vector<double> YErr_Data;

        // avg and stdv
        for(int i = 0; i < RAWyData.at(0).size(); i++) {
            std::vector<double> Y_temp_data;
            for (std::vector<double> Y_data_raw : RAWyData) {
                Y_temp_data.emplace_back(Y_data_raw.at(i));
            }
            std::tuple<double, double> mean_se = get_mean_and_se(Y_temp_data);
            Y_Data.emplace_back(std::get<0>(mean_se));
            YErr_Data.emplace_back(std::get<1>(mean_se));
        }
        hlp::saveOutData(filename + ".txt", header, x_lbl, y_lbl, xData, Y_Data, YErr_Data, N);
#ifdef SAVE_WITH_STEP_SIZE
        hlp::saveOutData("step_size_data/" + filename + "_step" + std::to_string(step) + ".txt",
                         header, x_lbl, y_lbl, xData, Y_Data, YErr_Data, N);
#endif

#endif
    }

    void saveAvgData(const std::string &filename, const std::string &header, const std::string &x_lbl, const std::string &y_lbl,
                         const std::vector<double> &xData, const std::vector<std::vector<double>> &RAWyData, const int &N) {

        std::vector<double> Y_Data;
        std::vector<double> YErr_Data;

        // avg and stdv
        for(int i = 0; i < RAWyData.at(0).size(); i++) {
            std::vector<double> Y_temp_data;
            for (std::vector<double> Y_data_raw : RAWyData) {
                Y_temp_data.emplace_back(Y_data_raw.at(i));
            }
            std::tuple<double, double> mean_se = get_mean_and_se(Y_temp_data);
            Y_Data.emplace_back(std::get<0>(mean_se));
            YErr_Data.emplace_back(std::get<1>(mean_se));
        }

        // save to file
        std::ofstream file;
        try {
            file.open(filename);
            file << header <<"\n\n";
            file << x_lbl << "\t" << y_lbl << "\t" << "yErr" << "\n";
            for (int i = 0; i < xData.size(); i++) {
                file << xData.at(i) << "\t" << Y_Data.at(i) << "\t" << YErr_Data.at(i) << "\n";
            }
        } catch (...) {
            file.close();
            std::cout << "failed to save to file\n";
        }
        file.close();
    }

}
