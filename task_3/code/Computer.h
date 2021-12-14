#ifndef COMPUTER_H
#define COMPUTER_H

#include <stdlib.h>
#include <cmath>
#include <vector>
#include <utility>

class Computer {

public:

    double a_1;
    double a_2;
    double b_1;
    double b_2;
    double h_1;
    double h_2;
    double alpha_t;
    double alpha_b;
    double alpha_l;
    double alpha_r;
    int64_t N;
    int64_t M;
    double** matrix_A;
    double* matrix_B;
    double* x_grid;
    double* y_grid;
    std::vector<std::pair<int64_t, int64_t> > index_map;

    double f_function(double x, double y);

    double q_function(double x, double y);

    double k_function(double x, double y);

    double left_border_function(double x, double y);

    double right_border_function(double x, double y);

    double top_border_function(double x, double y);

    double bottom_border_function(double x, double y);

    int64_t get_row_index(int64_t i, int64_t j);

    void create_matrices();

    double* fill_grid_vector(double first_element, double last_element, double step, int64_t count);

    Computer();

    Computer(double a_1, double b_1, double a_2, double b_2, double h_1, double h_2);

    Computer(double a_1, double b_1, double a_2, double b_2, int64_t N, int64_t M);

    double dot_product(double* u, double* v);

    double norm(double* u);

    void fill_index_map();

    double** from_vec_to_mat(double* vec);

    double** get_A();

    double* get_B();

    int64_t get_x_size();

    int64_t get_y_size();

    int64_t get_size();
    
    ~Computer();
};

#endif