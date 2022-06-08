#pragma once

#include <limits>

#include "../solver/base-solver.h"
#include "amr-node.h"
#include "amr-parameters.h"
#include "grid-boundary.h"
#include "refinement.h"

template <class SolverType, class ProblemType, class Fields, unsigned padding>
class AMRSolver : public Solver<AMRSolver<SolverType, ProblemType, Fields, padding>, Fields,
                                ProblemType, padding> {
  public:
    AMRSolver(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem,
              const AMRParameters& configuration);

    void init();
    void singleStep();
    void adjustConfig();

  private:
    const AMRParameters& configuration;
    Refinery<SolverType, ProblemType, Fields, padding> refinery;
    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>>& nodes;

    void levelStep(const unsigned level);
    void synchronise(const unsigned level);
    double minTimeDelta(const unsigned level);
    void updateTimeDelta(const double timeDelta, const unsigned level);
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRSolver<SolverType, ProblemType, Fields, padding>::AMRSolver(PaddedGrid<Fields, padding>& grid,
                                                               const Problem<ProblemType>& problem,
                                                               const AMRParameters& configuration)
    : Solver<AMRSolver<SolverType, ProblemType, Fields, padding>, Fields, ProblemType, padding>(
          grid, problem),
      configuration(configuration),
      refinery(Refinery<SolverType, ProblemType, Fields, padding>(problem, configuration)),
      nodes(refinery.getNodes()) {}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::init() {
    this->refinery.initialRefine(this->grid, this->problem);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::singleStep() {
    this->levelStep(0);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::synchronise(const unsigned level) {
    for (AMRNode<SolverType, ProblemType, Fields, padding>& node1 : this->nodes[level]) {
        for (AMRNode<SolverType, ProblemType, Fields, padding>* node2 : node1.siblings) {
            // Only want to exchange values once
            if (node1.id < node2->id) {
                std::optional<GridBoundary::Boundary> overlappingArea =
                    GridBoundary::overlappingArea(node1.gridBoundary, node2->gridBoundary);
                if (overlappingArea.has_value()) { // Should always be the case but you never know
                    GridBoundary::Boundary& overlap = overlappingArea.value();
                    for (unsigned x = overlap[0].first; x <= overlap[0].second; x++) {
                        const unsigned x1 = padding + x - node1.gridBoundary[0].first;
                        const unsigned x2 = padding + x - node2->gridBoundary[0].first;
                        for (unsigned y = overlap[1].first; y <= overlap[1].second; y++) {
                            const unsigned y1 = padding + y - node1.gridBoundary[1].first;
                            const unsigned y2 = padding + y - node2->gridBoundary[1].first;
                            for (unsigned z = overlap[2].first; z <= overlap[2].second; z++) {
                                const unsigned z1 = padding + z - node1.gridBoundary[2].first;
                                const unsigned z2 = padding + z - node2->gridBoundary[2].first;
                                for (unsigned field = 0; field < Fields().size(); field++) {
                                    node1.grid(x1, y1, z1)[field] =
                                        0.5 * node1.grid(x1, y1, z1)[field] +
                                        0.5 * node2->grid(x2, y2, z2)[field];
                                    node2->grid(x2, y2, z2)[field] = node1.grid(x1, y1, z1)[field];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::levelStep(const unsigned level) {
    if (level == 0) {
        for (AMRNode<SolverType, ProblemType, Fields, padding>& root : this->nodes[level]) {
            root.solver.step();
        }
        this->levelStep(level + 1);
        this->synchronise(level);
    } else if (level < this->nodes.size()) {
        for (unsigned step = 0; step < this->configuration.refinementFactor; step++) {
            for (AMRNode<SolverType, ProblemType, Fields, padding>& node : this->nodes[level]) {
                node.solver.step();
            }
            this->levelStep(level + 1);
            this->synchronise(level);
            if (step + 1 < this->configuration.refinementFactor) {
                for (AMRNode<SolverType, ProblemType, Fields, padding>& node : this->nodes[level]) {
                    node.updateBoundarySiblings();
                }
            }
        }
        for (AMRNode<SolverType, ProblemType, Fields, padding>& node : this->nodes[level]) {
            node.injectParents();
            node.updateBoundary();
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
double AMRSolver<SolverType, ProblemType, Fields, padding>::minTimeDelta(const unsigned level) {
    if (level < this->nodes.size()) {
        double minTimeDelta = this->minTimeDelta(level + 1) /
                              static_cast<double>(this->configuration.refinementFactor);
        for (AMRNode<SolverType, ProblemType, Fields, padding>& node : this->nodes[level]) {
            node.solver.adjust();
            minTimeDelta = std::min(minTimeDelta, node.solver.timeDelta);
        }
        return minTimeDelta;
    } else {
        return std::numeric_limits<double>::max();
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::updateTimeDelta(const double timeDelta,
                                                                          const unsigned level) {
    if (level < this->nodes.size()) {
        for (AMRNode<SolverType, ProblemType, Fields, padding>& node : this->nodes[level]) {
            node.solver.timeDelta = timeDelta;
        }
        updateTimeDelta(timeDelta / static_cast<double>(this->configuration.refinementFactor),
                        level + 1);
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::adjustConfig() {
    const double minTimeDelta = this->minTimeDelta(0);
    const double newTimeDelta = std::min(this->timeEnd - this->timeCurrent, minTimeDelta);
    this->timeDelta = newTimeDelta;
    this->updateTimeDelta(newTimeDelta, 0);
    if (this->timeStep + 1 % this->configuration.refinementInterval == 0) {
        this->refinery.refine();
    }
}
