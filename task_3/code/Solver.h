#ifndef SOLVER_H
#define SOLVER_H

#include <cstdlib>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "Computer.h"

const char FILENAME = "res.csv";

class Solver {
    int64_t omp_count;
    int64_t node_count;
    int64_t current_rank;
    int64_t grid_size;
    int64_t thread_count;
    std::vector<std::vector<int64_t>> grid;
    Computer computer;
    double EPS;
    double step;
    double* current_w;
    double* next_w;
    double* error_vec;
    double* multiplied_error;

public:

    Solver(Computer &computer, double eps);

    void solve();

    void write_result();

    ~Solver();
};

#endif