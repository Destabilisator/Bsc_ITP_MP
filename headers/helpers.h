#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <complex>
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

// invert first N bits of the integer s and returns result
int invertBits(int s, int N);

// sum up first N bits in s and returns result
int bitSum(int s, int N);

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

// provides representative of s to r and number of needed translations in l
void representative(int s, int *r, int *l, int N);

// provides representative of s to r and number of needed translations in l
// a possible reflection will be indicated by q
void representative(int s, int *r, int *l, int *q, int N);

/////////////////////////////// saving data ///////////////////////////////

// saves real (double) eigenvalues in ist eiVals to file filename and adds the header inscription above
void saveEiVals(const std::string &filename, const std::string &header, const std::list<double> &eiVals, const int &N);

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

// saves a vector of tuples (double) to be processed by python
void saveOutData(const std::string &filename, const std::string &header, const std::string &x_label,
                 const std::string &y_label, const std::vector<std::tuple<double, double>> &outData, const int &N);

/////////////////////////////// calculate quantities ///////////////////////////////

// returns the Matrix S^2 for a system of a given size
Eigen::MatrixXd spinMatrix(const int &N, const  int &SIZE);

// returns the Matrix S^2 for a given set of states
Eigen::MatrixXd spinMatrix(const int &N, const std::vector<int> &states);

// calculates the specific heat of a system for a set og eigenvalues and a beta
double getSpecificHeat(const double &beta, const std::vector<std::complex<double>>& eiVals, const int &N);

// calculates the magnetization of a system for a set og eigenvalues and a beta
double getSusceptibility(const double &beta, const Eigen::MatrixXcd &M, const std::vector<std::complex<double>>& eiVals, const int &N);

/////////////////////////////// others ///////////////////////////////

// check validity of user input
// [executable] N J_START J_END J_COUNT CORES SILENT
void validateInput(int argc, char* argv[], int *N, int *SIZE, double *J_START, double *J_END, int *J_COUNT,
                   const unsigned int *cpu_cnt, bool *silent, int *cores, const double *J1, const double *J2);
