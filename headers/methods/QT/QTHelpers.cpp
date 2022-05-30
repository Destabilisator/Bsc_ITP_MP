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
        std::cout << "saving to file '" << std::to_string(N) << "_" << filename << "'..." << std::endl;
        std::ofstream file;
        try {
            file.open("./results/" + std::to_string(N) + "_" + filename);
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

}