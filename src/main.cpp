#include "grid/padded-grid.h"
#include "parameters/constants.h"
#include "parameters/problem.h"
#include "solver/runge-kutta-solver.h"

int main(int argc, char** argv) {

    Problem problem = Problem::initialiseTestProblem();

    PaddedGrid<FieldStruct, GHOST_CELLS> grid({}, problem.numberCells[Direction::DirX],
                                              problem.numberCells[Direction::DirY],
                                              problem.numberCells[Direction::DirZ]);

    RungeKuttaSolver solver(grid, problem);

    std::cout << "------ Solving Grid ------" << std::endl << std::endl;
    solver.solve();
    std::cout << grid(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS)[0];

    return EXIT_SUCCESS;
}
