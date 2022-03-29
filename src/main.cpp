#include <cstdlib>
#include <iostream>
#include <memory>

#include "grid/padded-grid.h"
#include "grid/simple-grid.h"
#include "solver/simple-solver.h"

constexpr size_t CELLS_PER_DIMENSION = 10;

int main(int argc, char** argv) {

    SimpleGrid<double> grid = SimpleGrid<double>::initLinear(
        CELLS_PER_DIMENSION, CELLS_PER_DIMENSION, CELLS_PER_DIMENSION);
    grid.setBorderConst(3.0);
    grid.print();

    std::cout << "------ Solving Grid ------" << std::endl << std::endl;
    const double time = 1.0;
    const unsigned timeSteps = 10000;
    SimpleSolver solver(grid, time / timeSteps, time);
    solver.solve();

    grid.print();

    return EXIT_SUCCESS;
}
