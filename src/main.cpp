#include <cstdlib>
#include <iostream>
#include <memory>

#include "field-wrapper/phys-fields.h"
#include "grid/padded-grid.h"
#include "misc/constants.h"
#include "solver/runge-kutta-solver.h"

constexpr size_t CELLS_PER_DIMENSION = 10;

int main(int argc, char** argv) {

    PaddedGrid<FieldStruct, GHOST_CELLS> grid({}, CELLS_PER_DIMENSION, CELLS_PER_DIMENSION,
                                              CELLS_PER_DIMENSION);

    const double timeStart = 0.0;
    const double timeEnd = 1.0;
    const double timeDelta = 0.001;
    const unsigned rungeKuttaSteps = 2;
    RungeKuttaSolver solver(grid, timeDelta, timeStart, timeEnd, rungeKuttaSteps);

    std::cout << "------ Solving Grid ------" << std::endl << std::endl;
    solver.solve();
    std::cout << grid(CELLS_PER_DIMENSION / 2, CELLS_PER_DIMENSION / 2, CELLS_PER_DIMENSION / 2)[0];

    return EXIT_SUCCESS;
}
