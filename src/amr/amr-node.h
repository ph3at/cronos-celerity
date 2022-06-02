#pragma once

#include <cmath>
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
    std::optional<Fields> valueAtDown(const double xPos, const double yPos,
                                      const double zPos) const;
    std::optional<Fields> valueAtUp(const double xGlobal, const double yGlobal,
                                    const double zGlobal) const;
    void removeChildren();
    void addParent(AMRNode<SolverType, ProblemType, Fields, padding>* parent);
    void addSibling(AMRNode<SolverType, ProblemType, Fields, padding>* sibling);
    void addChild(AMRNode<SolverType, ProblemType, Fields, padding>* child);

    PaddedGrid<FieldStruct, GHOST_CELLS>& grid;
    std::array<std::pair<unsigned, unsigned>, 3> gridBoundary;

    Solver<SolverType, Fields, ProblemType, padding> solver;

  private:
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> parents;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> siblings;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> children;

    const AMRParameters& configuration;

    void inject();
    Fields interpolateDown(const double xPos, const double yPos, const double zPos) const;
    Fields interpolateUp(const double xLocal, const double yLocal, const double zLocal) const;
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>::AMRNode(PaddedGrid<Fields, padding>& grid,
                                                           const Problem<ProblemType>& problem,
                                                           const AMRParameters& configuration,
                                                           const double timeDelta,
                                                           const double timeCurrent)
    : grid(grid), solver(SolverType(grid, problem)), configuration(configuration) {
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
    for (AMRNode<SolverType, ProblemType, Fields, padding>* childGrid : this->children) {
        childGrid->integrate(this->solver.timeCurrent + this->solver.timeDelta);
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
                for (AMRNode<SolverType, ProblemType, Fields, padding>* child : this->children) {
                    const std::optional<Fields> childValues = child->valueAtDown(xPos, yPos, zPos);
                    if (childValues.has_value()) {
                        for (unsigned field = 0; field < Fields().size(); field++) {
                            this->grid(x, y, z)[field] = childValues.value()[field];
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
AMRNode<SolverType, ProblemType, Fields, padding>::valueAtDown(const double xPos, const double yPos,
                                                               const double zPos) const {
    if (xPos < this->grid.posLeft[Direction::DirX] || xPos > this->grid.posRight[Direction::DirX] ||
        yPos < this->grid.posLeft[Direction::DirY] || yPos > this->grid.posRight[Direction::DirY] ||
        zPos < this->grid.posLeft[Direction::DirZ] || zPos > this->grid.posRight[Direction::DirZ]) {
        return std::nullopt;
    } else {
        return std::make_optional(this->interpolateDown(xPos, yPos, zPos));
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Fields AMRNode<SolverType, ProblemType, Fields, padding>::interpolateDown(const double xPos,
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
    return fields;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
std::optional<Fields> AMRNode<SolverType, ProblemType, Fields, padding>::valueAtUp(
    const double xGlobal, const double yGlobal, const double zGlobal) const {
    const double xLocal = xGlobal - static_cast<double>(this->gridBoundary[0].first - padding);
    const double yLocal = yGlobal - static_cast<double>(this->gridBoundary[1].first - padding);
    const double zLocal = zGlobal - static_cast<double>(this->gridBoundary[2].first - padding);
    if (xLocal < 0.0 || xLocal >= static_cast<double>(this->grid.xDim()) || yLocal < 0.0 ||
        yLocal >= static_cast<double>(this->grid.yDim()) || zLocal < 0.0 ||
        zLocal >= static_cast<double>(this->grid.zDim())) {
        return std::nullopt;
    } else {
        const double shift = 1.0 / static_cast<double>(this->configuration.refinementFactor);
        return std::make_optional(this->interpolateUp(
            std::min(static_cast<double>(this->grid.xDim() - 1), std::max(0.0, xLocal - shift)),
            std::min(static_cast<double>(this->grid.yDim() - 1), std::max(0.0, yLocal - shift)),
            std::min(static_cast<double>(this->grid.zDim() - 1), std::max(0.0, zLocal - shift))));
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Fields AMRNode<SolverType, ProblemType, Fields, padding>::interpolateUp(const double xLocal,
                                                                        const double yLocal,
                                                                        const double zLocal) const {
    const double xLowD = std::floor(xLocal);
    const double xWeight = xLocal - xLowD;
    const unsigned xLow = static_cast<unsigned>(xLowD);
    const unsigned xHigh = static_cast<unsigned>(std::ceil(xLocal));
    const double yLowD = std::floor(yLocal);
    const double yWeight = yLocal - yLowD;
    const unsigned yLow = static_cast<unsigned>(yLowD);
    const unsigned yHigh = static_cast<unsigned>(std::ceil(yLocal));
    const double zLowD = std::floor(zLocal);
    const double zWeight = zLocal - zLowD;
    const unsigned zLow = static_cast<unsigned>(zLowD);
    const unsigned zHigh = static_cast<unsigned>(std::ceil(zLocal));
    Fields fields = { 0.0 };
    for (unsigned field = 0; field < fields.size(); field++) {
        const double xInterLowLow = xWeight * this->grid(xLow, yLow, zLow)[field] +
                                    (1.0 - xWeight) * this->grid(xHigh, yLow, zLow)[field];
        const double xInterHighLow = xWeight * this->grid(xLow, yHigh, zLow)[field] +
                                     (1.0 - xWeight) * this->grid(xHigh, yHigh, zLow)[field];
        const double xInterLowHigh = xWeight * this->grid(xLow, yLow, zHigh)[field] +
                                     (1.0 - xWeight) * this->grid(xHigh, yLow, zHigh)[field];
        const double xInterHighHigh = xWeight * this->grid(xLow, yHigh, zHigh)[field] +
                                      (1.0 - xWeight) * this->grid(xHigh, yHigh, zHigh)[field];
        const double yInterLow = yWeight * xInterLowLow + (1.0 - yWeight) * xInterHighLow;
        const double yInterHigh = yWeight * xInterLowHigh + (1.0 - yWeight) * xInterHighHigh;
        fields[field] = zWeight * yInterLow + (1.0 - zWeight) * yInterHigh;
    }
    return fields;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::removeChildren() {
    this->children.clear();
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::addParent(
    AMRNode<SolverType, ProblemType, Fields, padding>* parent) {
    this->parents.push_back(parent);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::addSibling(
    AMRNode<SolverType, ProblemType, Fields, padding>* sibling) {
    this->siblings.push_back(sibling);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::addChild(
    AMRNode<SolverType, ProblemType, Fields, padding>* child) {
    this->children.push_back(child);
}
