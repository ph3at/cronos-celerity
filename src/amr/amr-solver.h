#pragma once

#include "../solver/base-solver.h"
#include "amr-node.h"
#include "amr-parameters.h"

template <class SolverType, class ProblemType, class Fields, unsigned padding>
class AMRSolver
    : Solver<AMRSolver<SolverType, ProblemType, Fields, padding>, Fields, ProblemType, padding> {
  public:
    AMRSolver(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem,
              const AMRParameters& configuration);

    void initialise();
    void singleStep();
    void adjust();

  private:
    const AMRParameters& configuration;
    AMRNode<SolverType, ProblemType, Fields, padding> root;
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRSolver<SolverType, ProblemType, Fields, padding>::AMRSolver(PaddedGrid<Fields, padding>& grid,
                                                               const Problem<ProblemType>& problem,
                                                               const AMRParameters& configuration)
    : Solver<AMRSolver<SolverType, ProblemType, Fields, padding>, Fields, ProblemType, padding>(
          grid, problem),
      configuration(configuration) {
    this->root(grid, problem, configuration, problem.timeDelta, 0.0);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::initialise() {
    // Refine
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::singleStep() {
    this->root.step();
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::adjust() {
    // Adjust time
    // Refine
}
