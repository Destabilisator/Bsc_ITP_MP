#include "QTHelpers.h"

namespace QT::hlp {

    /////////////////////////////// Runge-Kutta 4 ///////////////////////////////

    Eigen::VectorXcd rungeKutta4Block(const Eigen::VectorXcd &vec, const Eigen::MatrixXcd &H, const double &step) {
        Eigen::VectorXcd k1 = - 0.5 * H * vec;
        Eigen::VectorXcd k2 = - 0.5 * H * (vec + step * 0.5 * k1);
        Eigen::VectorXcd k3 = - 0.5 * H * (vec + step * 0.5 * k2);
        Eigen::VectorXcd k4 = - 0.5 * H * (vec + step * k3);

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
