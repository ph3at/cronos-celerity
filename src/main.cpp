#include "configuration/constants.h"
#include "configuration/shock-tube.h"
#include "grid/grid-functions.h"
#include "grid/padded-grid.h"
#include "solver/runge-kutta-solver.h"

int main(int argc, char** argv) {

    ShockTube problem = ShockTube::initialiseTestProblem();

    PaddedGrid<FieldStruct, GHOST_CELLS> grid({}, problem.numberCells[Direction::DirX],
                                              problem.numberCells[Direction::DirY],
                                              problem.numberCells[Direction::DirZ]);
    problem.initialiseGrid(grid);

    RungeKuttaSolver<ShockTube> solver(grid, problem);

    std::cout << "------ Solving Grid ------" << std::endl << std::endl;
    solver.solve();
    std::cout << grid(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS)[0];

    return EXIT_SUCCESS;
}
