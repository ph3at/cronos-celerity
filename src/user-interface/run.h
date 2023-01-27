#pragma once

#include "../amr/amr-parameters.h"
#include "../amr/amr-solver.h"
#include "../configuration/constants.h"
#include "../configuration/shock-tube.h"
#include "../grid/grid-functions.h"
#include "../grid/padded-grid.h"
#include "../solver/runge-kutta-solver.h"
#include "../solver/runge-kutta-sycl-solver.h"

template <class ProblemType> void runAMR(const ProblemType& problem, const AMRParameters& amrConfig) {
    std::cout << std::endl
              << "----------------- Solving Grid using Adaptive Mesh Refinement -----------------" << std::endl
              << std::endl;

    const bool doOutput = true;

    AMRSolver<RungeKuttaSolver<ProblemType, FieldStruct, GHOST_CELLS>, ShockTube, FieldStruct, GHOST_CELLS> solver(
        problem, amrConfig, doOutput);
    solver.initialise();

    solver.solve();
    solver.report();
}

template <class ProblemType> void runRKS(const ProblemType& problem) {
    std::cout << std::endl
              << "----------------- Solving Grid using Runge-Kutta-Solver -----------------" << std::endl
              << std::endl;

    const unsigned rungeKuttaSteps = 2;
    const bool doOutput = true;
    RungeKuttaSolver<ProblemType, FieldStruct, GHOST_CELLS> solver(problem, rungeKuttaSteps, doOutput);
    solver.initialise();

    solver.solve();
    solver.report();
}

template <class ProblemType> void runSyclRKS(const ProblemType& problem) {
    std::cout << std::endl
              << "----------------- Solving Grid using Runge-Kutta-SYCL-Solver -----------------" << std::endl
              << std::endl;

    const unsigned rungeKuttaSteps = 2;
    RungeKuttaSyclSolver<ProblemType, FieldStruct, GHOST_CELLS> solver(problem, rungeKuttaSteps);
    solver.initialise();

    int timeSteps = 0;
    while (!solver.isFinished()) {
        std::cout << "--- Starting time step " << timeSteps << " ---" << std::endl;
        solver.step();
        solver.adjust();
        ++timeSteps;
    }
    solver.finalise();
    solver.report(timeSteps);
}
