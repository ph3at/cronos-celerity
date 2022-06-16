#pragma once

#include "../amr/amr-parameters.h"
#include "../amr/amr-solver.h"
#include "../configuration/constants.h"
#include "../configuration/shock-tube.h"
#include "../grid/grid-functions.h"
#include "../grid/padded-grid.h"
#include "../solver/runge-kutta-solver.h"

template <class ProblemType> void runAMR() {
    std::pair<PaddedGrid<FieldStruct, GHOST_CELLS>, ProblemType> problem =
        ShockTube::initialiseTestProblem();

    const bool doOutput = true;
    const AMRParameters amrConfig = { .refinementFactor = 4,
                                      .refinementInterval = 4,
                                      .bufferSize = 4,
                                      .efficiencyThreshold = 0.6,
                                      .truncationErrorTreshold = 1e-4 };
    AMRSolver<RungeKuttaSolver<ProblemType, FieldStruct, GHOST_CELLS>, ShockTube, FieldStruct,
              GHOST_CELLS>
        solver(problem.first, problem.second, amrConfig, doOutput);
    solver.initialise();
    std::cout << "----------------- Solving Grid -----------------" << std::endl << std::endl;

    solver.solve();
    solver.report();
}

template <class ProblemType> void runRK() {
    std::pair<PaddedGrid<FieldStruct, GHOST_CELLS>, ProblemType> problem =
        ShockTube::initialiseTestProblem();

    const unsigned rungeKuttaSteps = 2;
    const bool doOutput = true;
    RungeKuttaSolver<ProblemType, FieldStruct, GHOST_CELLS> solver(problem.first, problem.second,
                                                                   rungeKuttaSteps, doOutput);
    solver.initialise();
    std::cout << "----------------- Solving Grid -----------------" << std::endl << std::endl;

    solver.solve();
    solver.report();
}
