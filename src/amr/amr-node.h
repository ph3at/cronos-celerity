#pragma once

#include <optional>
#include <vector>

#include "../configuration/problem.h"
#include "../data-types/phys-fields.h"
#include "../solver/base-solver.h"
#include "amr-parameters.h"

template <class SolverType, class ProblemType, class Fields, unsigned padding> class AMRNode {
  public:
    AMRNode(PaddedGrid<Fields, padding> grid, const Problem<ProblemType>& problem,
            const AMRParameters& configuration,
            AMRNode<SolverType, ProblemType, Fields, padding>& parent);

    void integrate(const double untilTime);
    std::optional<Fields&> valueAt(const double xPos, const double yPos, const double zPos) const;

  private:
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>&> parents;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>&> children;

    PaddedGrid<FieldStruct, GHOST_CELLS> grid;
    const AMRParameters& configuration;
    Solver<SolverType, Fields, ProblemType, padding> solver;

    void inject();
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>::AMRNode(
    PaddedGrid<Fields, padding> grid, const Problem<ProblemType>& problem,
    const AMRParameters& configuration, AMRNode<SolverType, ProblemType, Fields, padding>& parent)
    : configuration(configuration) {
    this->grid = grid;
    this->parents();
    this->parents.push_back(parent);
    this->children();
    solver = SolverType(grid, problem);
    solver.timeDelta = parent.timeDelta / configuration.refinementFactor;
    solver.timeCurrent = parent.timeCurrent;
    solver.timeStep = 0;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::integrate(const double untilTime) {
    this->solver.timeEnd = untilTime;
    while (!this->solver.isFinished()) {
        for (AMRNode<SolverType, ProblemType, Fields, padding> childGrid : this->children) {
            childGrid.integrate(this->solver.timeCurrent + this->solver.timeDelta);
        }
        this->solver.step();
        this->inject();
    }
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
                double inChildren = 0.0;
                Fields childrenValues = { 0.0 };
                for (AMRNode<SolverType, ProblemType, Fields, padding>& child : this->children) {
                    const std::optional<Fields&> childValues = child.valueAt(xPos, yPos, zPos);
                    if (childValues.has_value()) {
                        inChildren += 1.0;
                        for (unsigned field = 0; field < Fields().size(); field++) {
                            childrenValues[field] += childValues.value()[field];
                        }
                    }
                }
                for (unsigned field = 0; field < Fields().size(); field++) {
                    this->grid(x, y, z)[field] =
                        0.5 * this->grid(x, y, z)[field] + 0.5 * childrenValues[field] / inChildren;
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::optional<Fields&>
AMRNode<SolverType, ProblemType, Fields, padding>::valueAt(const double xPos, const double yPos,
                                                           const double zPos) const {
    if (xPos < this->grid.posLeft[Direction::DirX] || xPos > this->grid.posRight[Direction::DirX] ||
        yPos < this->grid.posLeft[Direction::DirY] || yPos > this->grid.posRight[Direction::DirY] ||
        zPos < this->grid.posLeft[Direction::DirZ] || zPos > this->grid.posRight[Direction::DirZ]) {
        return std::nullopt;
    } else {
        // No interpolation for now
        const unsigned x =
            std::max(this->grid.xStart(),
                     std::min(this->grid.xEnd(),
                              static_cast<unsigned>((xPos - this->grid.posLeft[Direction::DirX]) *
                                                        grid.inverseCellSize[Direction::DirX] -
                                                    0.5) +
                                  padding));
        const unsigned y =
            std::max(this->grid.yStart(),
                     std::min(this->grid.yEnd(),
                              static_cast<unsigned>((yPos - this->grid.posLeft[Direction::DirY]) *
                                                        grid.inverseCellSize[Direction::DirY] -
                                                    0.5) +
                                  padding));
        const unsigned z =
            std::max(this->grid.zStart(),
                     std::min(this->grid.zEnd(),
                              static_cast<unsigned>((zPos - this->grid.posLeft[Direction::DirZ]) *
                                                        grid.inverseCellSize[Direction::DirZ] -
                                                    0.5) +
                                  padding));
        return std::make_optional(this->grid(x, y, z));
    }
}
