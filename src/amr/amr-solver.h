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
    AMRSolver(PaddedGrid<Fields, padding>& grid,
              const Problem<ProblemType, Fields, padding>& problem,
              const AMRParameters& configuration, const bool doOutput = false);

    void init();
    void singleStep();
    void adjustConfig();
    void finaliseResult();

  private:
    const AMRParameters& configuration;
    Refinery<SolverType, ProblemType, Fields, padding> refinery;
    std::vector<std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>>& nodes;

    void levelStep(const unsigned level);
    void synchronise(const unsigned level);
    double minTimeDelta(const unsigned level);
    void updateTimeDelta(const double timeDelta, const unsigned level);
    void extractGrid();
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRSolver<SolverType, ProblemType, Fields, padding>::AMRSolver(
    PaddedGrid<Fields, padding>& grid, const Problem<ProblemType, Fields, padding>& problem,
    const AMRParameters& configuration, const bool doOutput)
    : Solver<AMRSolver<SolverType, ProblemType, Fields, padding>, Fields, ProblemType, padding>(
          grid, problem, doOutput),
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
    double minTimeDelta = this->minTimeDelta(0);
    if (this->preciseEnd) {
        minTimeDelta = std::min(this->timeEnd - this->timeCurrent, minTimeDelta);
    }
    this->timeDelta = minTimeDelta;
    this->updateTimeDelta(minTimeDelta, 0);
    if ((this->timeStep + 1) % this->configuration.refinementInterval == 0 && !this->isFinished()) {
        this->refinery.refine();
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::finaliseResult() {
    for (std::vector<AMRNode<SolverType, ProblemType, Fields, padding>>& level : this->nodes) {
        for (AMRNode<SolverType, ProblemType, Fields, padding>& node : level) {
            node.solver.finalise();
        }
    }
    this->extractGrid();
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRSolver<SolverType, ProblemType, Fields, padding>::extractGrid() {
    const size_t sizeFactor = pow(this->configuration.refinementFactor, this->nodes.size() - 1);
    const size_t xSize = (this->grid.xEnd() - this->grid.xStart()) * sizeFactor;
    const size_t ySize = (this->grid.yEnd() - this->grid.yStart()) * sizeFactor;
    const size_t zSize = (this->grid.zEnd() - this->grid.zStart()) * sizeFactor;

    PaddedGrid<Fields, padding> result(this->grid.defaultValue, xSize, ySize, zSize,
                                       this->grid.posLeft, this->grid.posRight,
                                       this->grid.boundaryTypes);
    const double invRefinementFactor =
        1.0 / static_cast<double>(this->configuration.refinementFactor);
    for (unsigned xGlobal = 0, x = result.xStart(); xGlobal < xSize; xGlobal++, x++) {
        for (unsigned yGlobal = 0, y = result.yStart(); yGlobal < ySize; yGlobal++, y++) {
            for (unsigned zGlobal = 0, z = result.zStart(); zGlobal < zSize; zGlobal++, z++) {
                bool foundCell = false;
                for (AMRNode<SolverType, ProblemType, Fields, padding>& node : this->nodes.back()) {
                    if (node.isInBounds(xGlobal, yGlobal, zGlobal)) {
                        result(x, y, z) = node.valueAt(xGlobal, yGlobal, zGlobal);
                        foundCell = true;
                        break;
                    }
                }
                if (!foundCell) {
                    double xLevel = static_cast<double>(xGlobal) + 0.5;
                    double yLevel = static_cast<double>(yGlobal) + 0.5;
                    double zLevel = static_cast<double>(zGlobal) + 0.5;
                    for (int level = static_cast<int>(this->nodes.size()) - 2;
                         !foundCell && level >= 0; level--) {
                        xLevel *= invRefinementFactor;
                        yLevel *= invRefinementFactor;
                        zLevel *= invRefinementFactor;
                        for (AMRNode<SolverType, ProblemType, Fields, padding>& node :
                             this->nodes[level]) {
                            std::optional<Fields> fields = node.valueAtUp(xLevel, yLevel, zLevel);
                            if (fields.has_value()) {
                                result(x, y, z) = fields.value();
                                foundCell = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    this->grid.swap(result);
}
