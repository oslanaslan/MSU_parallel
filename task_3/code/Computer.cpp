#include "Computer.h"

double Computer::f_function(double x, double y) {
    // F(x, y) got from u(x, y), k(x, y), q(x, y)
    // using equation: F(x, y) = -\Delta u(x, y) + q(x, y)u(x, y)
    return 32 * x * y / std::pow(x * x + y * y + 1, 3) + 
        2 / (x * x + y * y + 1) +
            8 / std::pow(x * x + y * y + 1, 2);
}

double Computer::q_function(double x, double y) {
    // q(x, y) from task
    return 1;
}

double Computer::k_function(double x, double y) {
    // k(x, y) from task
    return 1 + std::pow(x + y, 2);
}

double Computer::left_border_function(double x, double y) {
    // x == A1
    return this->f_function(this->a_1, y);
}

double Computer::right_border_function(double x, double y) {
    // x == A2
    return this->f_function(this->a_2, y);
}

double Computer::top_border_function(double x, double y) {
    // y == B2
    return this->f_function(x, this->b_2);
}

double Computer::bottom_border_function(double x, double y) {
    // y == B1
    return this->f_function(x, this->b_1);
}

int64_t Computer::get_row_index(int64_t i, int64_t j) {
    // Flattern vector indices
    return (this->M + 1) * i + j;
}

void Computer::create_matrices() {
    // Fill A and B matrices for solving Aw = B equation
    // where w is flattern matrix of w_ij grid elements in vector form
    if (this->matrix_A && this->matrix_B) {
        return;
    }

    int64_t a_size = (this->N + 1) * (this->M + 1);
    this->matrix_A = new double*[a_size];
    this->matrix_B = new double[a_size];

    for (int64_t i = 0; i < a_size; i++) {
        this->matrix_A[i] = new double[a_size];
        this->matrix_B[i] = 0;

        for (int64_t j = 0; j < a_size; j++) {
            this->matrix_A[i][j] = 0;
        }
    }

    for (int64_t i = 0; i < this->N; i++) {
        for (int64_t j = 0; j < this->M; j++) {
            int64_t row_i = this->get_row_index(i, j);
            double x = this->x_grid[i];
            double x_plus = x + this->h_1 / 2;
            double x_minus = x - this->h_1 / 2;
            double y = this->y_grid[j];
            double y_plus = y + this->h_2 / 2;
            double y_minus = y - this->h_2 / 2;
            double h_1_mult = 1 / (this->h_1 * this->h_1);
            double h_2_mult = 1 / (this->h_2 * this->h_2);

            if (i != 0 && i != this->M && j != 0 && j != this->N) {
                // Inner points of grid
                
                double xi_1 = h_1_mult * this->k_function(x_plus, y) +
                    h_1_mult * this->k_function(x_minus, y) +
                        h_2_mult * this->k_function(x, y_plus) +
                            h_2_mult * this->k_function(x, y_minus) +
                                this->q_function(x, y);
                double xi_2 = h_1_mult * this->k_function(x_plus, y);
                double xi_3 = h_1_mult * this->k_function(x_minus, y);
                double xi_4 = h_2_mult * this->k_function(x, y_plus);
                double xi_5 = h_2_mult * this->k_function(x, y_minus);
                
                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i + 1, j)] = -xi_2;
                this->matrix_A[row_i][this->get_row_index(i - 1, j)] = -xi_3;
                this->matrix_A[row_i][this->get_row_index(i, j + 1)] = -xi_4;
                this->matrix_A[row_i][this->get_row_index(i, j - 1)] = -xi_5;

                this->matrix_B[row_i] = this->f_function(x, y);
            }
            else if (i == 0 && j != 0 && j != this->N) {
                // Left border points
                x_minus = x_grid[i + 1] - this->h_1 / 2;
                double xi_1 = this->q_function(x, y) + 2 * this->alpha_l / this->h_1 +
                    2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1) - this->k_function(x, y_minus) / this->h_2;
                double xi_2 = -2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1);
                double xi_3 = this->k_function(x, y_minus) / this->h_2;
                
                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i + 1, j)] = xi_2;
                this->matrix_A[row_i][this->get_row_index(i, j - 1)] = xi_3;

                this->matrix_B[row_i] = this->f_function(i, j) + 2 * this->left_border_function(x, y) / this->h_1;
            }
            else if (i == this->M && j != 0 && j != this->N) {
                // Right border points
                double xi_1 = this->q_function(x, y) + 2 * this->alpha_r / this->h_1 +
                    2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1) - this->k_function(x, y_minus) / this->h_2;
                double xi_2 = this->k_function(x, y_minus) / this->h_2;
                double xi_3 = 2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1);

                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i, j - 1)] = xi_2;
                this->matrix_A[row_i][this->get_row_index(i - 1, j)] = -xi_3;

                this->matrix_B[row_i] = this->f_function(x, y) + 2 * this->right_border_function(x, y) / this->h_1;
            }
            else if (i != 0 && i != this->N && j == 0) {
                // Bottom border points
                y_minus = this->y_grid[j + 1] - this->h_2 / 2;

                double xi_1 = 2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2) + this->q_function(x, y) +
                    2 * this->alpha_b / this->h_2 - this->k_function(x_minus, y) / this->h_1;
                double xi_2 = this->k_function(x_minus, y) / this->h_1;
                double xi_3 = -2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2);

                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i - 1, j)] = xi_2;
                this->matrix_A[row_i][this->get_row_index(i, j + 1)] = xi_3;

                this->matrix_B[row_i] = this->f_function(x, y) + 2 * this->top_border_function(x, y) / this->h_2;

            }
            else if (i != 0 && i != this->M && j == this->N) {
                // Top border points
                double xi_1 = 2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2) + this->q_function(x, y) +
                    2 * this->alpha_t / this->h_2 - this->k_function(x_minus, y) / this->h_1;
                double xi_2 = -2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2);
                double xi_3 = this->k_function(x_minus, y) / this->h_1;

                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i, j - 1)] = xi_2;
                this->matrix_A[row_i][this->get_row_index(i - 1, j)] = xi_3;

                this->matrix_B[row_i] = this->f_function(x, y) + 2 * this->top_border_function(x, y) / this->h_2;
            }
            else if (i == 0 && j == 0) {
                // P(0, 0) corner
                x_minus = this->x_grid[i + 1] - this->h_1 / 2;
                y_minus = this->y_grid[j + 1] - this->h_2 / 2;

                double xi_1 = 2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1) +
                    2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2) + this->q_function(x, y) + 
                        2 * this->alpha_l / this->h_1 + 2 * this->alpha_b / this->h_2;
                double xi_2 = -2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1);
                double xi_3 = -2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2);

                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i + 1, j)] = xi_2;
                this->matrix_A[row_i][this->get_row_index(i, j + 1)] = xi_3;

                this->matrix_B[row_i] = this->f_function(x, y) + (2 / this->h_1 + 2 / this->h_2) * this->left_border_function(x, y);
                
            }
            else if (i == 0 && j == this->N) {
                // P(0, N) corner
                x_minus = this->x_grid[i + 1] - this->h_1 / 2;

                double xi_1 = 2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1) +
                    2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2) + this->q_function(x, y) + 
                        2 * this->alpha_l / this->h_1 + 2 * this->alpha_t / this->h_2;
                double xi_2 = -2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1);
                double xi_3 = -2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2);

                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i + 1, j)] = xi_2;
                this->matrix_A[row_i][this->get_row_index(i, j - 1)] = xi_3;

                this->matrix_B[row_i] = this->f_function(x, y) + (2 / this->h_1 + 2 / this->h_2) * this->right_border_function(x, y);
            }
            else if (i == this->M && j == 0) {
                // P(M, 0) corner
                y_minus = this->y_grid[j + 1] - this->h_2 / 2;

                double xi_1 = 2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1) +
                    2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2) + this->q_function(x, y) + 
                        2 * this->alpha_r / this->h_1 + 2 * this->alpha_b / this->h_2;
                double xi_2 = -2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1);
                double xi_3 = -2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2);

                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i - 1, j)] = xi_2;
                this->matrix_A[row_i][this->get_row_index(i, j + 1)] = xi_3;

                this->matrix_B[row_i] = this->f_function(x, y) + (2 / this->h_1 + 2 / this->h_2) * this->right_border_function(x, y);
            }
            else if (i == this->M && j == this->N) {
                // P(M, N) corner
                double xi_1 = 2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1) +
                    2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2) + this->q_function(x, y) + 
                        2 * this->alpha_r / this->h_1 + 2 * this->alpha_t / this->h_2;
                double xi_2 = -2 * this->k_function(x_minus, y) / (this->h_1 * this->h_1);
                double xi_3 = -2 * this->k_function(x, y_minus) / (this->h_2 * this->h_2);

                this->matrix_A[row_i][row_i] = xi_1;
                this->matrix_A[row_i][this->get_row_index(i + 1, j)] = xi_2;
                this->matrix_A[row_i][this->get_row_index(i, j - 1)] = xi_3;

                this->matrix_B[row_i] = this->f_function(x, y) + (2 / this->h_1 + 2 / this->h_2) * this->right_border_function(x, y);
            }
        }
    }
}

double* Computer::fill_grid_vector(double first_element, double last_element, double step, int64_t count) {
    // Fill given array of length count with elements from first to last with step
    double* vector = new double[count];

    vector[0] = first_element;
    for (int64_t i = 1; i < count; i++) {
        vector[i] = vector[i - 1] + step;
    }
    vector[count - 1] = last_element;

    return vector;
}

Computer::Computer() {
    this->a_1 = -2;
    this->a_2 = 3;
    this->b_1 = -1;
    this->b_2 = 4;
    this->h_1 = 0.01;
    this->h_2 = 0.01;
    // TODO check if this returns correct number of grid elements
    this->N = (a_2 - a_1) / h_1;
    this->M = (b_2 - b_1) / h_2;
    this->matrix_A = NULL;
    this->matrix_B = NULL;
    this->alpha_l = 0;
    this->alpha_r = 0;
    this->alpha_t = 1;
    this->alpha_b = 1;
    //TODO check if generates correct grid
    this->x_grid = this->fill_grid_vector(a_1, a_2, h_1, this->N + 1);
    this->y_grid = this->fill_grid_vector(b_1, b_2, h_2, this->M + 1);
    this->create_matrices();
    this->fill_index_map();
}

Computer::Computer(double a_1, double b_1, double a_2, double b_2, double h_1, double h_2) {
    this->a_1 = a_1;
    this->a_2 = a_2;
    this->b_1 = b_1;
    this->b_2 = b_2;
    this->h_1 = h_1;
    this->h_2 = h_2;
    // TODO check if this returns correct number of grid elements
    this->N = (a_2 - a_1) / h_1;
    this->M = (b_2 - b_1) / h_2;
    this->matrix_A = NULL;
    this->matrix_B = NULL;
    this->alpha_l = 0;
    this->alpha_r = 0;
    this->alpha_t = 1;
    this->alpha_b = 1;
    //TODO check if generates correct grid
    this->x_grid = this->fill_grid_vector(a_1, a_2, h_1, this->N + 1);
    this->y_grid = this->fill_grid_vector(b_1, b_2, h_2, this->M + 1);
    this->create_matrices();
    this->fill_index_map();
}

Computer::Computer(double a_1, double b_1, double a_2, double b_2, int64_t N, int64_t M) {
    this->a_1 = a_1;
    this->a_2 = a_2;
    this->b_1 = b_1;
    this->b_2 = b_2;
    // TODO check this
    this->h_1 = (a_2 - a_1) / N;
    this->h_2 = (b_2 - b_1) / M;
    this->N = N;
    this->M = M;
    this->matrix_A = NULL;
    this->matrix_B = NULL;
    this->alpha_l = 0;
    this->alpha_r = 0;
    this->alpha_t = 1;
    this->alpha_b = 1;
    //TODO check if generates correct grid
    this->x_grid = this->fill_grid_vector(a_1, a_2, h_1, this->N + 1);
    this->y_grid = this->fill_grid_vector(b_1, b_2, h_2, this->M + 1);
    this->create_matrices();
    this->fill_index_map();
}

double Computer::dot_product(double* u, double* v) {
    // Count dot special dot product in H space of linear grid functions
    double res = 0;

    for (int64_t i = 0; i < this->N; i++ ) {
        double sub_res = 0;

        for (int64_t j = 0; j < this->M; j++) {
            double weight = (i > 0 && i < N ? 1 : 0.5) * (j > 0 && j < M ? 1 : 0.5);
            sub_res += weight * this->h_2 * u[this->get_row_index(i, j)] * v[this->get_row_index(i, j)];
        }

        res += this->h_1 * sub_res;
    }

    return res;
}

double Computer::norm(double* u) {
    // Norm in space with given dot product
    return sqrt(this->dot_product(u, u));
}

void Computer::fill_index_map() {
    for (int64_t i = 0; i < this->N + 1; i++) {
        for (int64_t j = 0; j < this->M + 1; j++) {
            this->index_map[this->get_row_index(i, j)] = std::make_pair(i, j);
        }
    }
}

double** Computer::from_vec_to_mat(double* vec) {
    // Convert W vec to W matrix
    double** mat = new double*[this->N + 1];

    for (int64_t i = 0; i < this->N + 1; i++) {
        mat[i] = new double[this->M + 1];

        for (int64_t j = 0; j < this->M + 1; j++) {
            mat[i][j] = vec[this->get_row_index(i, j)];
        }
    }

    return mat;
}

double** Computer::get_A() {
    return this->matrix_A;
}

double* Computer::get_B() {
    return this->matrix_B;
}

int64_t Computer::get_x_size() {
    return this->N + 1;
}

int64_t Computer::get_y_size() {
    return this->M + 1;
}

int64_t Computer::get_size() {
    return (this->N + 1) * (this->M + 1);
}

Computer::~Computer() {
    if (!this->matrix_A) {
        int64_t a_size = (this->N + 1) * (this->M + 1);

        for (int64_t i = 0; i < a_size; i++) {
            delete[] this->matrix_A[i];
        }

        delete[] this->matrix_A;
    }
    if (!this->matrix_B) {
        delete[] this->matrix_B;
    }
    if (!this->x_grid) {
        delete[] this->x_grid;
    }
    if (!this->y_grid) {
        delete[] this->y_grid;
    }
}