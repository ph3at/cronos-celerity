#pragma once

#include <optional>
#include <vector>

#include "../configuration/problem.h"
#include "../data-types/phys-fields.h"
#include "../solver/base-solver.h"
#include "amr-parameters.h"

template <class SolverType, class ProblemType, class Fields, unsigned padding> class AMRNode {
  public:
    AMRNode(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem,
            const AMRParameters& configuration, const double timeDelta, const double timeCurrent);

    void integrate(const double untilTime);
    void step();
    std::optional<Fields> valueAt(const double xPos, const double yPos, const double zPos) const;

  private:
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>&> parents;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>&> children;

    PaddedGrid<FieldStruct, GHOST_CELLS>& grid;
    const AMRParameters& configuration;
    Solver<SolverType, Fields, ProblemType, padding> solver;

    void inject();
    Fields interpolate(const double xPos, const double yPos, const double zPos) const;
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>::AMRNode(PaddedGrid<Fields, padding>& grid,
                                                           const Problem<ProblemType>& problem,
                                                           const AMRParameters& configuration,
                                                           const double timeDelta,
                                                           const double timeCurrent)
    : configuration(configuration) {
    this->grid = grid;
    this->parents();
    this->children();
    this->solver = SolverType(grid, problem);
    this->solver.timeDelta = timeDelta;
    this->solver.timeCurrent = timeCurrent;
    this->solver.timeStep = 0;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::integrate(const double untilTime) {
    this->solver.timeEnd = untilTime;
    while (!this->solver.isFinished()) {
        this->step();
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::step() {
    for (AMRNode<SolverType, ProblemType, Fields, padding> childGrid : this->children) {
        childGrid.integrate(this->solver.timeCurrent + this->solver.timeDelta);
    }
    this->solver.step();
    this->inject();
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::inject() {
    double xPos = this->grid.posLeft[Direction::DirX];
    double yPos = this->grid.posLeft[Direction::DirX];
    double zPos = this->grid.posLeft[Direction::DirX];
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd();
         x++, xPos += this->grid.cellSize[Direction::DirX]) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd();
             y++, yPos += this->grid.cellSize[Direction::DirY]) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd();
                 z++, zPos += this->grid.cellSize[Direction::DirZ]) {
                // Could/Should be optimised to try last child with value first.
                // Potentially switch to child calling parent
                for (AMRNode<SolverType, ProblemType, Fields, padding>& child : this->children) {
                    const std::optional<Fields> childValues = child.valueAt(xPos, yPos, zPos);
                    if (childValues.has_value()) {
                        for (unsigned field = 0; field < Fields().size(); field++) {
                            this->grid(x, y, z)[field] = childValues.value[field];
                        }
                        break;
                    }
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::optional<Fields>
AMRNode<SolverType, ProblemType, Fields, padding>::valueAt(const double xPos, const double yPos,
                                                           const double zPos) const {
    if (xPos < this->grid.posLeft[Direction::DirX] || xPos > this->grid.posRight[Direction::DirX] ||
        yPos < this->grid.posLeft[Direction::DirY] || yPos > this->grid.posRight[Direction::DirY] ||
        zPos < this->grid.posLeft[Direction::DirZ] || zPos > this->grid.posRight[Direction::DirZ]) {
        return std::nullopt;
    } else {
        return std::make_optional(this->interpolate(xPos, yPos, zPos));
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Fields AMRNode<SolverType, ProblemType, Fields, padding>::interpolate(const double xPos,
                                                                      const double yPos,
                                                                      const double zPos) const {
    const unsigned xStart = static_cast<unsigned>((xPos - this->grid.posLeft[Direction::DirX]) *
                                                  grid.inverseCellSize[Direction::DirX]) +
                            padding;
    const unsigned yStart = static_cast<unsigned>((yPos - this->grid.posLeft[Direction::DirY]) *
                                                  grid.inverseCellSize[Direction::DirY]) +
                            padding;
    const unsigned zStart = static_cast<unsigned>((zPos - this->grid.posLeft[Direction::DirZ]) *
                                                  grid.inverseCellSize[Direction::DirZ]) +
                            padding;
    Fields fields = { 0.0 };
    for (unsigned x = xStart; x < xStart + this->configuration.refinementFactor; x++) {
        for (unsigned y = yStart; y < yStart + this->configuration.refinementFactor; y++) {
            for (unsigned z = zStart; z < zStart + this->configuration.refinementFactor; z++) {
                for (unsigned field = 0; field < fields.size(); field++) {
                    fields[field] += this->grid(x, y, z)[field];
                }
            }
        }
    }
    const double invSqrRefinementFactor =
        1.0 / (static_cast<double>(this->configuration.refinementFactor) *
               static_cast<double>(this->configuration.refinementFactor));
    for (unsigned field = 0; field < fields.size(); field++) {
        fields[field] *= invSqrRefinementFactor;
    }
}
