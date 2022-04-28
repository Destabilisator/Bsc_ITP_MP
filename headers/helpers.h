#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <complex>

/////////////////////////////// bits ///////////////////////////////

// print out first N bits of the integer s
void printBits(int s, int N);

// cyclicly translate first N bits in s by n to the right and returns result
int translateRight(int s, int n, int N);

// cyclicly translate first N bits in s by n to the left and returns result
int translateLeft(int s, int n, int N);

// reflect first N bits of the integer s and returns result
int reflectBits(int s, int N);

// invert first N bits of the integer s and returns result
int invertBits(int s, int N);

// sum up first N bits in s and returns result
int bitSum(int s, int N);

/////////////////////////////// states ///////////////////////////////

// finds all possible states with magnetization m_z = m and outputs to states
void fillStates(std::vector<int> *states, int m, int N, int size);

// returns the index of the state s in the vector states
int findState(std::vector<int> *states, int s, int M);

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
void saveEiVals(const std::string &filename, const std::string &header, const std::list<double> &eiVals);

// saves complex (std::complex<double>) eigenvalues in ist eiVals to file filename and adds the header inscription above
void saveComplexEiVals(const std::string &filename, const std::string &header, const std::list<std::complex<double>> &eiVals);

// saves a real matrix (double**) of size size wo the file filename and adds the header inscription above
void saveHamilton(double** hamilton, const std::string &filename, const std::string &header, int size);

// saves a complex matrix (std::complex<double>**) of size size wo the file filename and adds the header inscription above
void saveComplexHamilton(std::complex<double> **hamilton,const std::string &filename, const std::string &header, int size);
