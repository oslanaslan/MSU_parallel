#include "Computer.h"
#include "Solver.h"

double EPS = 1e-5;

int main() {
    Computer computer;
    Solver solver(computer, EPS);

    solver.solve();
    solver.write_result();

    return 0;
}