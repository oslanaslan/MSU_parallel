#include Solver.h

Solver::Solver(Computer &computer, double eps) {
    // MPI_Init();
    // MPI_Comm_size(MPI_COMM_WORLD, &this->node_count);
    // MPI_Comm_rank(MPI_COMM_WORLD, &this->current_rank);
    this->thread_count = omp_get_max_threads();
    this->EPS = eps;
    this->step = 0;
    this->computer = computer;
    this->grid_size = computer.get_size();
    this->grid.resize(this->node_count);
    this->current_w = new double[this->grid_size];
    this->next_w = new double[this->grid_size];
    this->error_vec = new double[this->grid_size];
    this->multiplied_error = new double[this->grid_size];

    int64_t min_batch_size = this->grid_size / this->thread_count;
    std::vector<int64_t> vector_sizes;

    vector_sizes.resize(this->thread_count);

    for (int64_t i = 0; i < this->thread_count; i++) {
        vector_sizes[i] = min_batch_size;
    }

    for (int64_t i = min_batch_size * this->thread_count; i < this->grid_size; i++) {
        vector_sizes[i - min_batch_size * this->thread_count] += 1;
    }

    int64_t index = 0;

    for (int64_t thread = 0; thread < this->thread_count; thread++) {
        this->grid[thread].resize(vector_sizes[thread]);

        for (int64_t i = 0; i < vector_sizes[thread]; i++) {
            this->grid[thread][i] = index++;
        }
    }

    std::fill(this->current_w, this->current_w + this->grid_size, 0);
    std::fill(this->next_w, this->next_w + this->grid_size, 0);
    std::fill(this->error_vec, this->error_vec + this->grid_size, 0);
}

void Solver::solve() {
    double tolerance = 100 * this->EPS;

    do
    {
        #pragma omp parallel for
        for (int64_t thread = 0; thread < this->thread_count; thread++) {
            for (int64_t index : this->grid[thread]) {
                this->next_w[index] = this->current_w[index] - this->step * this->error_vec[index];
            }

            #pragma omp barrier

            for (int64_t index : this->grid[thread]) {
                this->error_vec[index] = 0;

                for (int64_t i = 0; i < this->grid_size; i++) {
                    this->error_vec[index] += this->computer.matrix_A[index][i] * this->next_w[i] - this->computer.matrix_B[i];
                }
            }

            #pragma omp barrier

            for (int64_t index : this->grid[thread]) {
                this->multiplied_error[index] = 0;

                for (int64_t i = 0; i < this->grid_size; i++) {
                    this->multiplied_error[index] += this->computer.matrix_A[index][i] * this->error_vec[i];
                }
            }

            #pragma omp barrier

            if (thread == 0) {
                // Main thread
                double dot_prod = this->computer.dot_product(this->multiplied_error, this->error_vec);
                double norm = this->computer.norm(this->multiplied_error);
                this->step = dot_prod / norm;
                tolerance = this->computer.norm(this->next_w, this->current_w);
            }

            #pragma omp barrier

            for (int64_t index : this->grid[thread]) {
                this-current_w[index] = this->next_w[index];
            }

            #pragma omp barrier
        }
    } while (tolerance > this->EPS);
}

void Solver::write_result() {
    double** mat = this->computer.from_vec_to_mat(this->current_w);
    std::ofstream file;
    file.open(FILENAME);

    for (int64_t i = 0; i < this->computer.N + 1; i++) {
        for (int64_t j = 0; j < this->computer.M + 1; j++) {
            file << mat[i][j] << ",";
        }
        file << "\n";
    }

    file.close();
}

Solver::~Solver() {
    delete[] current_w;
    delete[] next_w;
    delete[] error_vec;
    delete[] multiplied_error;
}