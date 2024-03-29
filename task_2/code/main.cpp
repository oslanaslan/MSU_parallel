#include "header.h"

using namespace std;

int main(int argc, char **argv) {
	int rank = -1;
	int iter_count = 1;
	bool is_tol_enough = false;
	double start_time = 0;
	double prog_time = 0;
	double count_value = 0;
	double tol = 0;
	double current_sum = 0;
	double res_value = 0;
	double all_sum = 0;
	double eps = atof(argv[1]);
	double max_prog_time = 0;

	// Init	
	MPI_Comm comm_gr;
	MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NODE_COUNT);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(rank + 1);
	
	// Start
	start_time = MPI_Wtime();

	while (!is_tol_enough && iter_count < MAX_ITER_COUNT) {
		// Generate sum of F func values
		current_sum = generate_values();

		MPI_Reduce(&current_sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

		if (rank == MASTER) {
			// Master
			count_value += all_sum;
			res_value = VOLUME * count_value / iter_count;
			tol = abs(I - res_value);

			if (tol < eps) {
				is_tol_enough = true;
			}
		}

		iter_count++;
	}
	
	prog_time = MPI_Wtime() - start_time;

	if (!is_tol_enough) {
		printf("Terminated due to max loop limit\n");
	}
	else {
		printf("Terminated due to good sollution\n");
	}

	// Count time
	MPI_Reduce(&prog_time, &max_prog_time, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);

	// Display result
	if (rank == MASTER) {
		printf("Proc count: %d\nEPS: %.10lf\n", NODE_COUNT, eps);
		printf("Res: %.10lf\nTol: %.10lf\nTime: %.10lf\n\n", res_value, tol, max_prog_time);
	}

	// Exit
	MPI_Finalize();

	return 0;
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
		x = (double)rand() / RAND_MAX;
		y = (double)rand() / RAND_MAX;
		z = (double)rand() / RAND_MAX;
		current_sum += f_func(x, y, z);
	}

	return current_sum / (NODE_COUNT * MAX_SAMPLES_COUNT);
}
