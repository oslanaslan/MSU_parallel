#include <cstdio>
#include <cstdlib>
#include "mpi.h"
#include <cmath>
#include <stdint.h>
#include <cstring>

#define MASTER		0

const double A_1 = 0;
const double A_2 = 0;
const double A_3 = 0;
const double B_1 = 1;
const double B_2 = 1;
const double B_3 = 1;
const double I = (double)log(2.0) / 2.0 - 5.0 / 16.0;
double VOLUME = (B_1 - A_1) * (B_2 - A_2) * (B_3 - A_3);
const uint64_t MAX_SAMPLES_COUNT = 10;
const uint64_t MAX_ITER_COUNT = 10000000;

int NODE_COUNT = 0;

double generate_values();
double f_func(double x, double y, double z);

