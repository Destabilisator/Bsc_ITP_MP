#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <complex>
#include <chrono>
#include <Eigen/Eigenvalues>

/////////////////////////////// bits ///////////////////////////////

// print out first N bits of the integer s
void printBits(int s, int N);

// cyclically translate first N bits in s by n to the right and returns result
int translateRight(int s, int n, int N);

// cyclically translate first N bits in s by n to the left and returns result
int translateLeft(int s, int n, int N);

// reflect first N bits of the integer s and returns result
int reflectBits(int s, int N);

// inverts first N bits of the integer s and returns result
int invertBits(int s, int N);

// sum up first N bits in s and returns result
int bitSum(int s, int N);

// inverts the first N bits of an integer s
//int invertBits(int s, int N);

/////////////////////////////// states ///////////////////////////////

// finds all possible states with magnetizationBlocksAnsatz m_z = m and outputs to states
void fillStates(std::vector<int> *states, int m, int N, int size);

// returns the index of the state s in the vector states
int findState(const std::vector<int>& states, int s);

int findState(const std::vector<std::tuple<int, int>>& states, int s);

// if s is the smallest state integer, returns its periodicity else -1
int checkState(int s, int k, int N);

// if s is the smallest state integer, writes periodicity in translation to r (without) and m (with) reflection
// if a smaller state exists r = -1
void checkState(int s, int *r, int *m, int k, int N);

// checks whether a smaller in can be reached through translation r, reflection mp or inversion mz
void checkStateSI(const int &s, int &r, int &mp, int &mz, int &mpz, const int &k, const int &N);

// provides representative of s to r and number of needed translations in l
void representative(int s, int *r, int *l, int N);

// provides representative of s to r and number of needed translations in l
// a possible reflection will be indicated by q
void representative(int s, int *r, int *l, int *q, int N);

// provides representative of s to r and number of needed translations in l
// a possible reflection will be indicated by q
// a possible spin inversion will be indicated by g
void representative(const int &s, int &r, int &l, int &q, int &g, const int &N);

/////////////////////////////// saving data ///////////////////////////////

// saves real (double) eigenvalues in ist eiVals to file filename and adds the header inscription above
void saveEiVals(const std::string &filename, const std::string &header, const std::list<double> &eiVals, const int &N);

void saveEiVals(const std::string &filename, const std::string &header, const std::vector<double> &eiVals, const int &N);

// saves complex (std::complex<double>) eigenvalues in ist eiVals to file filename and adds the header inscription above
void saveComplexEiVals(const std::string &filename, const std::string &header, const std::list<std::complex<double>> &eiVals, const int &N);

void saveComplexEiVals(const std::string &filename, const std::string &header, const std::vector<std::complex<double>> &eiVals, const int &N);

// saves a real matrix (double**) of size to the file filename and adds the header inscription above
void saveHamilton(double** hamilton, const std::string &filename, const std::string &header, const int &size, const int &N);

// saves a complex matrix (std::complex<double>**) of size size to the file filename and adds the header inscription above
void saveComplexHamilton(std::complex<double> **hamilton,const std::string &filename, const std::string &header, const int &size, const int &N);

// saves a real matrix (Eigen::MatrixXd) of arbitrary size to the file filename and adds the header inscription above
void saveMatrixToFile(const Eigen::MatrixXd& matrix, const std::string &filename, const std::string &header, const int &N);

// saves a complex matrix (Eigen::MatrixXcd) of arbitrary size to the file filename and adds the header inscription above
void saveComplexMatrixToFile(const Eigen::MatrixXcd& matrix, const std::string &filename, const std::string &header, const int &N);

// saves a vector of tuples (double, double) to be processed by python
void saveOutData(const std::string &filename, const std::string &header, const std::string &x_label,
                 const std::string &y_label, const std::vector<std::tuple<double, double>> &outData, const int &N);

// saves a vector of tuples (int, double) to be processed by python
void saveOutData(const std::string &filename, const std::string &header, const std::string &x_label,
                 const std::string &y_label, const std::vector<std::tuple<int, double>> &outData, const int &N);

// saves a vector of tuples (double, double) to be processed into a 3D-plot by python
void save3DPlotDataC(const double &J, const int &N, const std::vector<std::tuple<double, double>>& C_func_T);

void save3DPlotDataX(const double &J, const int &N, const std::vector<std::tuple<double, double>>& X_func_T);

/////////////////////////////// calculate quantities ///////////////////////////////

// returns the Matrix S^2 for a system of a given size
Eigen::MatrixXd spinMatrix(const int &N, const  int &SIZE);

// returns the Matrix S^2 for a given set of states
Eigen::MatrixXd spinMatrix(const int &N, const std::vector<int> &states);

// calculates the specific heat of a system for a set og eigenvalues and a beta
double getSpecificHeat(const double &beta, const std::vector<std::complex<double>>& eiVals, const int &N);

double getSpecificHeat(const double &temp, const std::vector<double>& eiVals, const int &N);

// calculates the magnetization of a system for a set og eigenvalues and a beta
double getSusceptibility(const double &beta, const Eigen::MatrixXcd &M, const std::vector<std::complex<double>>& eiVals, const int &N);

double getSusceptibilityDegeneracy(const double &temp, const Eigen::MatrixXcd &M, const std::vector<std::complex<double>>& eiVals, const int &N);

double getSusceptibilityDegeneracy(const double &temp, const Eigen::MatrixXd &M, const std::vector<double>& eiVals, const int &N);

double getSusceptibilityDegeneracy(const double &temp, const std::vector<Eigen::MatrixXd> &M_list, const std::vector<std::vector<double>>& eiVal_list, const int &N);

/////////////////////////////// others ///////////////////////////////

// check validity of user input
// [executable] N J_START J_END J_COUNT CORES SILENT
void validateInput(int &argc, char* argv[], const unsigned int &cpu_cnt, int &N, int &SIZE, double &J_START, double &J_END,
                   int &J_COUNT, double &T_START, double &T_END, int &T_COUNT, bool &silent, int &cores, bool &plotsIn3D,
                   bool skipSilent, const double &J1, const double &J2);

// converts time in s in hh:mm:ss
std::string formatTime(std::chrono::duration<double> elapsed_seconds);
