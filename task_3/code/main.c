#include "header.h"

int main(int argc, char **argv) {
	int rank = -1;
	int size = -1;
	int N = -1;
	int *data_vec = NULL;
	int data_vec_size = 0;
	int source = -1;

	// Init	
	MPI_Comm comm_gr;
	MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Read data
	data_vec = read_data(argc, argv, rank, size, &data_vec_size, &source);
	if (data_vec == NULL) { 
		return 0;
	}

	// Share data
	N = floor(log2(size));
	int start_ind = 0;
	for (; rank >= round(pow(2, start_ind)); start_ind++); 
	int src_num = N - start_ind;
	int addrs[src_num];
	int inp_addrs[src_num];
	data_vec = share_data(data_vec, &data_vec_size, rank, size, &start_ind, &src_num, addrs, inp_addrs, N);
	
	// Sort
	qsort_alg(data_vec, 0, data_vec_size-1);
	
	// Receive sorted data
	int *sorted_addrs[src_num];
	int inp_sizes[src_num];
	data_vec = recv_data(data_vec, &data_vec_size, rank, sorted_addrs, inp_sizes, addrs, inp_addrs, &src_num);

	// Save result
	if (rank == 0) {
		FILE *f = fopen("output.txt", "w");
		for (int i = 0; i < data_vec_size; i++)
			fprintf(f, "%d ", data_vec[i]);
		fprintf(f, "\n");	
		fclose(f);
		for (int i = 0; i < data_vec_size; i++)
			printf("%d ", data_vec[i]);
		printf("\n");
		my_bcast(NO_DATA, size);
	}
	else {
		MPI_Send(&data_vec_size, 1, MPI_INT, source, DATA_TAG, MPI_COMM_WORLD);
		MPI_Send(data_vec, data_vec_size, MPI_INT, source, DATA_TAG, MPI_COMM_WORLD);
	}

	// Exit
	free(data_vec);
	MPI_Finalize();
	return 0;
}

void qsort_alg(int *numbers, int left, int right) {
	int pivot; 
	int l_hold = left; 
	int r_hold = right; 

	pivot = numbers[left];
	while (left < right) {
		while ((numbers[right] >= pivot) && (left < right))
			right--; 
		if (left != right) {
			numbers[left] = numbers[right]; 
			left++; 
		}
		while ((numbers[left] <= pivot) && (left < right))
			left++; 
		if (left != right) {
			numbers[right] = numbers[left]; 
			right--; 
		}
	}
	numbers[left] = pivot; 
	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot) 
		qsort_alg(numbers, left, pivot - 1);
	if (right > pivot)
		qsort_alg(numbers, pivot + 1, right);
}

void my_bcast(int msg, int size) {
	for (int i = 0; i < size; i++)
		MPI_Send(&msg, 1, MPI_INT, i, MSG_TAG, MPI_COMM_WORLD);
}


int *read_data(int argc, char **argv, int rank, int size, int *data_vec_size_ret, int *source_ret) {
	int data_vec_size = 0;
	int *data_vec = NULL;
	int source = -1;
	MPI_Status stat;

	if (rank == 0) {
			// rank 0
			if (argc == 2) {
				int file_len = 0;
				int mum;
				FILE *input = fopen(argv[1], "rt");
				if (input != NULL) {
					int num;	
					while (fscanf(input, "%d", &num) == 1) {
						data_vec_size++;
						data_vec = (int *)realloc(data_vec, sizeof(int)*data_vec_size);
						data_vec[data_vec_size-1] = num;
					}
					fclose(input);
					my_bcast(OK_MSG, size);
				}
				else {
					printf("No such file: %s\n", argv[1]);
					my_bcast(ERR_MSG, size);
					free(data_vec);
					MPI_Finalize();
					return NULL;
				}
			}
			else if (argc > 2) {
				int num;
				data_vec_size = argc - 1;
				data_vec = (int *)realloc(data_vec, sizeof(int)*data_vec_size);
				for (int i = 1; i < argc; i++) {
					data_vec[i-1] = atoi(argv[i]);
				}
				
				my_bcast(OK_MSG, size);
			}
			else {
				printf("No data\n");
				my_bcast(ERR_MSG, size);
				free(data_vec);
				MPI_Finalize();	
				return NULL;
			}
		}
		else {
			// rank 1..N-1
			int my_stat = -1;
			MPI_Recv(&my_stat, 1, MPI_INT, ROOT, MSG_TAG, MPI_COMM_WORLD, &stat);
			if (my_stat == ERR_MSG) {
				free(data_vec);
				MPI_Finalize();
				return NULL;
			}	
			MPI_Recv(&data_vec_size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
			if (data_vec_size == NO_DATA) {
				free(data_vec);
				MPI_Finalize();
				return NULL;
			}
			data_vec = (int *)malloc(sizeof(int)*data_vec_size);
			MPI_Recv(data_vec, data_vec_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
			source = stat.MPI_SOURCE;

		}
	
	*data_vec_size_ret = data_vec_size;	
	*source_ret = source;
	return data_vec;
}

int *share_data(int *data_vec, int *data_vec_size_ret, int rank, int size, int *start_ind_ret, int *src_num_ret, int *addrs, int *inp_addrs, int N) {
	int data_vec_size = *data_vec_size_ret;
	int start_ind = *start_ind_ret;
	int src_num = *src_num_ret;

	
	// Share data
	// 0 -> 0 + 2^0, 0 + 2^1, ... , 0 + 2^(N-1)
	// 1 -> 1 + 2^1, 1 + 2^2, ... , 1 + 2^(N-1)
	// 2 -> 2 + 2^2, 2 + 2^3, ... , 2 + 2^(N-1)
	// ...

	for (int i = start_ind; i < N; i++) {
		addrs[i - start_ind] = rank + round(pow(2, i));		
		int addr = rank + round(pow(2, i));
		int marker = data_vec[data_vec_size / 2];
		int less_vec_size = 0;
		int more_vec_size = 0;
		int *less_vec = NULL;
		int *more_vec = NULL;
		for (int j = 0; j < data_vec_size; j++) {
			if (data_vec[j] < marker) {
				less_vec = (int *)realloc(less_vec, sizeof(int) * (++less_vec_size));
				less_vec[less_vec_size-1] = data_vec[j];
			}
			else {
				more_vec = (int *)realloc(more_vec, sizeof(int) * (++more_vec_size));
				more_vec[more_vec_size-1] = data_vec[j];
			}
		}
		if (more_vec_size == 0 || less_vec_size == 0) {
			src_num = i - start_ind;
			free(less_vec);
			free(more_vec);
			for (int j = i; j < N; j++) {
				int addr = rank + round(pow(2, j));
				int msg = NO_DATA;
				MPI_Send(&msg, 1, MPI_INT, addr, DATA_TAG, MPI_COMM_WORLD);
			}
			break;
		}
		else {
		MPI_Send(&more_vec_size, 1, MPI_INT, addr, DATA_TAG, MPI_COMM_WORLD);		
		MPI_Send(more_vec, more_vec_size, MPI_INT, addr, DATA_TAG, MPI_COMM_WORLD);
		free(data_vec);
		free(more_vec);
		data_vec = less_vec;
		data_vec_size = less_vec_size;
		}
	}
	
	*data_vec_size_ret = data_vec_size;
	*start_ind_ret = start_ind;
	*src_num_ret = src_num;
	return data_vec;
}

int *recv_data(int *data_vec, int *data_vec_size_ret, int rank, int **sorted_addrs_ret, int *inp_sizes_ret, int *addrs_ret, int *inp_addrs_ret, int *src_num_ret) {
	int **sorted_addrs = sorted_addrs_ret;
	int *inp_sizes = inp_sizes_ret;
	int data_vec_size = *data_vec_size_ret;
	int *addrs = addrs_ret;
	int *inp_addrs = inp_addrs_ret;
	int src_num = *src_num_ret;
	MPI_Status stat;

	for (int i = 0; i < src_num; i++) {
		int vec_size = 0;
		MPI_Recv(&vec_size, 1, MPI_INT, MPI_ANY_SOURCE, DATA_TAG, MPI_COMM_WORLD, &stat);
		sorted_addrs[i] = (int *)malloc(sizeof(int) * vec_size);
		MPI_Recv(sorted_addrs[i], vec_size, MPI_INT, stat.MPI_SOURCE, DATA_TAG, MPI_COMM_WORLD, &stat);
		inp_addrs[i] = stat.MPI_SOURCE;
		inp_sizes[i] = vec_size;
	}	
	for (int i = src_num-1; i >= 0; i--) {
		int j = 0;
		while (addrs[i] != inp_addrs[j]) 
			j++;
		data_vec = (int *)realloc(data_vec, sizeof(int) * (data_vec_size + inp_sizes[j]));
		int *bigger_data_vec = sorted_addrs[j];
		for (int k = 0; k < inp_sizes[j]; k++)
			data_vec[k + data_vec_size] = bigger_data_vec[k];
		data_vec_size += inp_sizes[j];
		free(sorted_addrs[j]);
	}

	*data_vec_size_ret = data_vec_size;
	*src_num_ret = src_num;
	return data_vec;
}
