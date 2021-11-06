#include <cstdio>
#include <cstdlib>
#include "mpi.h"
#include <cmath>
#include <time.h>
#include <stdint.h>

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
const double I = log(2) / 2 - 5 / 16;
double VOLUME = (B_1 - A_1) * (B_2 - A_2) * (B_3 - A_3);
const double EPS = 1e-5;
const char RESULT_FILEPATH[100] = "results.txt";
const uint64_t MAX_SAMPLES_COUNT = 1000;
int NODE_COUNT = 10;
const uint64_t MAX_ITER_COUNT = 1000000;

struct RecResult {
    int type;
	double res;
	bool success;
	uint64_t time;
    int msg;
};

bool send_msg(int source, int msg);
bool send_data(int source, double result, uint64_t cur_time);
void save_data(double* tol_vec, uint64_t max_time, double* res_vec, uint64_t iter_count);
double generate_values();
double f_func(double x, double y, double z);
RecResult recv_msg(int source);
RecResult recv_data(int source);
RecResult process_data(double current_sum, double* tol_vec, double* res_vec, RecResult slaves_res, uint64_t current_iter);



