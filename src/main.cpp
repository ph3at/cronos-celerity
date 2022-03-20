#include <iostream>
#include <stdlib.h>

#include "grid/simple-grid.h"
#include "solver/simple-solver.h"

using namespace std;

const size_t CELLS_PER_DIMENSION = 10;

void solve(Grid* grid);

int main(int argc, char** argv) {

    SimpleGrid grid =
        SimpleGrid::init_linear(CELLS_PER_DIMENSION, CELLS_PER_DIMENSION, CELLS_PER_DIMENSION);
    grid.set_border_const(3.0);
    grid.print();

    cout << "------ Solving Grid ------" << endl << endl;
    solve(&grid);
    grid.print();

    return EXIT_SUCCESS;
}

void solve(Grid* grid) {
    SimpleSolver simple_solver = SimpleSolver();
    Solver* solver = &simple_solver;
    const unsigned time_steps = 10000;
    const double time_step = 1.0 / double(time_steps);
    for (unsigned step = 0; step < time_steps; step++) {
        solver->compute_step(grid, time_step);
    }
}
