#include <cstdio>
#include <cstdlib>
#include "mpi.h"
#include <cmath>
#include <time.h>
#include <stdint.h>
#include <cstring>

// typedef int bool;
// enum {false, true};

#define ERR_MSG		1
#define OK_MSG		2
#define CONTINUE    3
#define EXIT        4
#define DATA_TAG	1
#define MSG_TAG		2
#define MASTER		0
#define ALL         -1
#define ANY         MPI_ANY_SOURCE
// Tag 1 - data, 2 - msg

const double A_1 = 0;
const double A_2 = 0;
const double A_3 = 0;
const double B_1 = 1;
const double B_2 = 1;
const double B_3 = 1;
const double I = (double)log(2.0) / 2.0 - 5.0 / 16.0;
double VOLUME = (B_1 - A_1) * (B_2 - A_2) * (B_3 - A_3);
double EPS = 1e-3;
const char RESULT_FILEPATH[100] = "results.txt";
const uint64_t MAX_SAMPLES_COUNT = 10000;
int NODE_COUNT = 10;
const uint64_t MAX_ITER_COUNT = 1000;

double generate_values();
double f_func(double x, double y, double z);
