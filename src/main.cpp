#include <cstdlib>
#include <iostream>
#include <memory>

#include "grid/simple-grid.h"
#include "solver/simple-solver.h"

constexpr size_t CELLS_PER_DIMENSION = 10;

void solve(Grid& grid);

int main(int argc, char** argv) {

    SimpleGrid grid =
        SimpleGrid::initLinear(CELLS_PER_DIMENSION, CELLS_PER_DIMENSION, CELLS_PER_DIMENSION);
    grid.setBorderConst(3.0);
    grid.print();

    std::cout << "------ Solving Grid ------" << std::endl << std::endl;
    solve(grid);
    grid.print();

    return EXIT_SUCCESS;
}

void solve(Grid& grid) {
    SimpleSolver simpleSolver = SimpleSolver();
    std::unique_ptr<Solver> solver = std::make_unique<SimpleSolver>(simpleSolver);
    const unsigned timeSteps = 10000;
    const double time = 1.0;
    const double deltaTime = time / double(timeSteps);
    for (unsigned step = 0; step < timeSteps; step++) {
        solver->computeStep(grid, deltaTime);
    }
}
