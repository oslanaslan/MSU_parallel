#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <time.h>

#define ERR_MSG		1
#define OK_MSG		2
#define DATA_TAG	1
#define MSG_TAG		2
#define ROOT		0
#define NO_DATA		-1
// Tag 1 - data, 2 - msg

void qsort_alg(int *numbers, int left, int right);
void my_bcast(int msg, int size);
int *read_data(int argc, char **argv, int rank, int size, int *data_vec_size_ret, int *source_ret);
int *share_data(int *data_vec, int *data_vec_size, int rank, int size, int *start_ind, int *src_num, int *addrs, int *inp_addrs, int N);
int *recv_data(int *data_vec, int *data_vec_size, int rank, int **sorted_addrs, int *inp_sizes, int *addrs, int *inp_addrs, int *src_num);
