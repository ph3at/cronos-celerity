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

    RungeKuttaSolver<ShockTube> solver(grid, problem);
    solver.initialise();
    std::cout << "----------------- Solving Grid -----------------" << std::endl << std::endl;
    solver.solve();
    solver.report();

    const SimpleGrid<FieldStruct> resultOfOld =
        GridFunctions::readFromFile("test-data/step-16.dat");
    GridFunctions::compare(resultOfOld, grid, true);

    return EXIT_SUCCESS;
}
