#include "header.h"

int main(int argc, char **argv) {
	int rank = -1;
	double tol_vec[MAX_ITER_COUNT];
	double res_vec[MAX_ITER_COUNT];
	int source = -1;
	bool is_tol_enough = false;
	time_t start_time = 0;
	uint64_t iter_count = 0;

	printf("%d: started\n", rank);

	// Init	
	MPI_Comm comm_gr;
	MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NODE_COUNT);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	start_time = time(NULL);
	srand(start_time * (rank + 1));

	printf("%d: init: success", rank);

	while (!is_tol_enough && iter_count < MAX_ITER_COUNT) {
		printf("%d: while: %ld iteration started\n", rank, iter_count);
		// Generate sum of F func values
		double current_sum = generate_values();

		if (rank == MASTER) {
			// Master
			struct RecResult res = recv_data(ALL);

			if (res.success) {
				struct RecResult res = process_data(current_sum, tol_vec, res_vec, res, iter_count);
				is_tol_enough = res.success;
			}
			if (is_tol_enough) {
				uint64_t max_time = (uint64_t)(time(NULL) - start_time);

				max_time = max_time > res.time ? max_time : res.time;
				save_data(tol_vec, max_time, res_vec, iter_count);
				send_msg(ALL, EXIT);
			}
			else {
				send_msg(ALL, CONTINUE);
			}
		}
		else {
			// Slave
			send_data(MASTER, current_sum, (uint64_t)(time(NULL) - start_time));
			struct RecResult res = recv_msg(MASTER);

			if (!res.success || res.type != MSG_TAG && res.msg == EXIT) {
				is_tol_enough = true;
			}
		}

		iter_count++;
	}

	// Exit
	if (rank == MASTER) {
		send_msg(ALL, EXIT);
	}
	MPI_Finalize();

	return 0;
}

bool send_msg(int source, int msg) {
	bool res = true;

	if (source != ALL) {
		res &= MPI_Send(&msg, 1, MPI_INT, source, MSG_TAG, MPI_COMM_WORLD) ? true : false;
	}
	else {
		for (int i = 0; i < NODE_COUNT; i++) {
			res &= MPI_Send(&msg, 1, MPI_INT, i, MSG_TAG, MPI_COMM_WORLD) ? true : false;
		}
	}

	return res;
}

bool send_data(int source, double result, uint64_t cur_time) {
	bool res = true;

	if (source != ALL) {
		res &= MPI_Send(&result, 1, MPI_DOUBLE, source, DATA_TAG, MPI_COMM_WORLD) ? true : false;
		res &= MPI_Send(&cur_time, 1, MPI_UINT64_T, source, DATA_TAG, MPI_COMM_WORLD) ? true : false;
	}
	else {
		for (int i = 0; i < NODE_COUNT; i++) {
			res &= MPI_Send(&result, 1, MPI_DOUBLE, i, DATA_TAG, MPI_COMM_WORLD) ? true : false;
			res &= MPI_Send(&cur_time, 1, MPI_UINT64_T, i, DATA_TAG, MPI_COMM_WORLD) ? true : false;
		}
	}

	return res;
}

struct RecResult recv_msg(int source) {
	struct RecResult result;
	MPI_Status stat;
	int msg;

	MPI_Recv(&msg, 1, MPI_INT, source, MSG_TAG, MPI_COMM_WORLD, &stat);

	result.type = MSG_TAG;
	result.msg = msg;
	result.success = true;

	return result;
}

struct RecResult recv_data(int source) {
	struct RecResult result;
	MPI_Status stat;
	double data;
	uint64_t cur_time;

	if (source == ALL) {
		for (int i = 1; i < NODE_COUNT; i++) {
			MPI_Recv(&data, 1, MPI_DOUBLE, source, DATA_TAG, MPI_COMM_WORLD, &stat);
			MPI_Recv(&cur_time, 1, MPI_UINT64_T, source, DATA_TAG, MPI_COMM_WORLD, &stat);

			result.type = DATA_TAG;
			result.res += data;
			result.success = true;
			result.time = cur_time > result.time ? cur_time : result.time;
		}
	}
	else {
		MPI_Recv(&data, 1, MPI_DOUBLE, source, DATA_TAG, MPI_COMM_WORLD, &stat);
		MPI_Recv(&cur_time, 1, MPI_UINT64_T, source, DATA_TAG, MPI_COMM_WORLD, &stat);

		result.type = DATA_TAG;
		result.res = data;
		result.time = cur_time;
		result.success = true;
	}

	return result;
}

void save_data(double* tol_vec, uint64_t max_time, double* res_vec, uint64_t iter_count) {
	FILE *input = fopen(RESULT_FILEPATH, "rt");

	if (input) {
		fprintf(input, "{");
		fprintf(input, "\"time\": %ld", max_time);
		fprintf(input, "\"results\":[");
		
		for (uint64_t i = 0; i < iter_count + 1; i++) {
			fprintf(input, "%f,", res_vec[i]);
		}

		fprintf(input, "],\"tolerance\":[");

		for (uint64_t i = 0; i < iter_count + 1; i++) {
			fprintf(input, "%f,", tol_vec[i]);
		}

		fprintf(input, "]}");
		fclose(input);
	}
}

struct RecResult process_data(double current_sum, double* tol_vec, double* res_vec, struct RecResult slaves_res, uint64_t current_iter) {
	struct RecResult result;
	double cur_tol;

	result.res = VOLUME * (current_sum + slaves_res.res) / (NODE_COUNT * MAX_SAMPLES_COUNT);
	result.time = slaves_res.time;
	result.type = DATA_TAG;
	cur_tol = (I - result.res);
	result.success = cur_tol < EPS ? true : false;

	tol_vec[current_iter] = cur_tol;
	res_vec[current_iter] = result.res;

	return result;
}

double f_func(double x, double y, double z) {
	double res = 0;

	if (x >= 0 && y >= 0 && z >= 0 && z <= 1 - x - y) {
		res = 1 + x + y + z;
		res = 1 / (res * res * res);
	}
	
	return res;
}

double generate_values() {
	double current_sum = 0;
	double x, y, z;

	for (uint64_t i = 0; i < MAX_SAMPLES_COUNT; i++) {
		x = (B_1 - A_1) * rand() / RAND_MAX + A_1;
		y = (B_2 - A_2) * rand() / RAND_MAX + A_2;
		z = (B_3 - A_3) * rand() / RAND_MAX + A_3;
		current_sum += f_func(x, y, z);
	}

	return current_sum;
}