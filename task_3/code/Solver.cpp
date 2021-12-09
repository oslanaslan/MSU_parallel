#include <stdlib.h>

class Solver {

    double h_1;
    double h_2;
    int N;
    int M;
    double** matrix_A;
    double* matrix_B;

    double f_function(double x, double y) {
        // TODO
        return 0;
    }

    double q_function(double x, double y) {
        // TODO
        return 0;
    }

    double k_function(double x, double y) {
        // TODO
        return 0;
    }

    double left_border_function(double x, double y) {
        // TODO
        return 0;
    }

    double right_border_function(double x, double y) {
        // TODO
        return 0;
    }

    double upper_border_function(double x, double y) {
        // TODO
        return 0;
    }

    double bottom_border_function(double x, double y) {
        // TODO
        return 0;
    }

    int64_t get_row_index_for_w(int64_t i, int64_t j) {
        return (this->N + 1) * i + j;
    }

    void create_matrices() {
        int64_t a_size = (this->N + 1) * (this->M + 1);
        this->matrix_A = new double*[a_size];
        this->matrix_B = new double[a_size];

        for (int64_t i = 0; i < a_size; i++) {
            this->matrix_A[i] = new double[a_size];

            for (int64_t j = 0; i < a_size; i++) {
                // A[i][j] = ...

            }

            // B[i] = ...

        }


    }

public:

    Solver(double h_1, double h_2, int N, int M) {
        this->h_1 = h_1;
        this->h_2 = h_2;
        this->N = N;
        this->M = M;
        this->matrix_A = nullptr;
        this->matrix_B = nullptr;
    }

    double* from_matrix_to_vec(double** target_matrix) {
        // TODO
        return nullptr;
    }

    double** from_vec_to_matrix(double* target_vec) {
        // TODO
        return nullptr;
    }
    
    ~Solver() {
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
    }
};