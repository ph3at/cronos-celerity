#include "configuration/constants.h"
#include "configuration/shock-tube.h"
#include "grid/grid-functions.h"
#include "grid/padded-grid.h"
#include "solver/runge-kutta-solver.h"

int main(int argc, char** argv) {

    std::pair<PaddedGrid<FieldStruct, GHOST_CELLS>, ShockTube> shockTube =
        ShockTube::initialiseTestProblem();

    RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS> solver(shockTube.first, shockTube.second);
    solver.initialise();
    std::cout << "----------------- Solving Grid -----------------" << std::endl << std::endl;

    for (unsigned timeStep = 1; timeStep <= 16; timeStep++) {
        solver.step();
        const SimpleGrid<FieldStruct> baseline =
            GridFunctions::readFromFile("test-data/step-" + std::to_string(timeStep) + ".dat");
        GridFunctions::compare(baseline, shockTube.first, true, false);
    }
    solver.report();

    return EXIT_SUCCESS;
}
