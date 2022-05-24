#pragma once

#include "../solver/base-solver.h"
#include "amr-node.h"
#include "amr-parameters.h"
#include "refinement.h"

template <class SolverType, class ProblemType, class Fields, unsigned padding>
class AMRSolver : public Solver<AMRSolver<SolverType, ProblemType, Fields, padding>, Fields,
                                ProblemType, padding> {
  public:
    AMRSolver(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem,
              const AMRParameters& configuration);

    void initialise() {}
    void singleStep();
    void adjust();

  private:
    const AMRParameters& configuration;
    Refinery<SolverType, ProblemType, Fields, padding> refinery;
    AMRNode<SolverType, ProblemType, Fields, padding>& root;
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRSolver<SolverType, ProblemType, Fields, padding>::AMRSolver(PaddedGrid<Fields, padding>& grid,
                                                               const Problem<ProblemType>& problem,
                                                               const AMRParameters& configuration)
    : Solver<AMRSolver<SolverType, ProblemType, Fields, padding>, Fields, ProblemType, padding>(
          grid, problem),
      configuration(configuration),
      refinery(Refinery<SolverType, ProblemType, Fields, padding>(configuration.refinementFactor)),
      root(this->refinery.initialRefine(grid, problem, configuration)) {
    // this->refinery(configuration.refinementFactor);
    // this->root = refinery.initialRefine(this->grid, this->problem, this->configuration);
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
