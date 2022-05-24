#include "amr/amr-parameters.h"
#include "amr/amr-solver.h"
#include "configuration/constants.h"
#include "configuration/shock-tube.h"
#include "grid/grid-functions.h"
#include "grid/padded-grid.h"
#include "solver/runge-kutta-solver.h"

int main(int argc, char** argv) {

    std::pair<PaddedGrid<FieldStruct, GHOST_CELLS>, ShockTube> shockTube =
        ShockTube::initialiseTestProblem();

    const AMRParameters amrConfig = { .refinementFactor = 5,
                                      .refinementInterval = 5,
                                      .bufferSize = 4,
                                      .efficiencyThreshold = 0.6,
                                      .truncationErrorTreshold = 0.0005 };
    AMRSolver<RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS>, ShockTube, FieldStruct,
              GHOST_CELLS>
        solver(shockTube.first, shockTube.second, amrConfig);
    solver.initialise();
    std::cout << "----------------- Solving Grid -----------------" << std::endl << std::endl;

    solver.solve();
    solver.report();

    return EXIT_SUCCESS;
}
