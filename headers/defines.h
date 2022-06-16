#define PI 3.14159265358979323846
#define EPSILON 1e-10
//#define fixedPrecision // remove if not debug

// nested multithreading loops
#define OUTERMOST_NESTED_THREADS 3
#define OUTER_NESTED_THREADS 1
#define INNER_NESTED_THREADS 1

// ED data calculation
#define calc_C_X_over_beta

// QT data output
#define SAVE_WITH_SETS_OF_n_SAMPLES
#define SAVE_WITH_DATA_FROM_ALL_SAMPLES
#define SAVE_WITH_STEP_SIZE

// prevent multithreading in Eigen
//#define EIGEN_DONT_PARALLELIZE
