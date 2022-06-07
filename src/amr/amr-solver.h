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
    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>>& nodes;

    void levelStep(const unsigned level);
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRSolver<SolverType, ProblemType, Fields, padding>::AMRSolver(PaddedGrid<Fields, padding>& grid,
                                                               const Problem<ProblemType>& problem,
                                                               const AMRParameters& configuration)
    : Solver<AMRSolver<SolverType, ProblemType, Fields, padding>, Fields, ProblemType, padding>(
          grid, problem),
      configuration(configuration),
      refinery(Refinery<SolverType, ProblemType, Fields, padding>(problem, configuration)),
      nodes(this->refinery.initialRefine(grid, problem)) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::singleStep() {
    this->levelStep(0);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::levelStep(const unsigned level) {
    if (level == 0) {
        for (AMRNode<SolverType, ProblemType, Fields, padding> root : this->nodes[level]) {
            root.solver.step();
        }
        this->levelStep(level + 1);
        // synchronise
    } else if (level < this->nodes.size()) {
        for (unsigned step = 0; step < this->configuration.refinementFactor; step++) {
            for (AMRNode<SolverType, ProblemType, Fields, padding> node : this->nodes[level]) {
                node.solver.step();
            }
            this->levelStep(level + 1);
        }
        // get boundaries/synchronise
        for (AMRNode<SolverType, ProblemType, Fields, padding> node : this->nodes[level]) {
            node.injectParents();
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::adjust() {
    // Adjust time
    this->refinery.refine();
}
