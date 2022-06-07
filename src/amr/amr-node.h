#pragma once

#include <cmath>
#include <optional>
#include <vector>

#include "../configuration/problem.h"
#include "../data-types/phys-fields.h"
#include "../solver/base-solver.h"
#include "amr-parameters.h"
#include "grid-boundary.h"

template <class SolverType, class ProblemType, class Fields, unsigned padding> class AMRNode {
  public:
    AMRNode(PaddedGrid<Fields, padding>& grid, const Problem<ProblemType>& problem,
            const AMRParameters& configuration, const double timeDelta, const double timeCurrent);

    void setValue(const unsigned xGlobal, const unsigned yGlobal, const unsigned zGlobal,
                  const Fields& fields);
    void injectParents();
    Fields valueAtDown(const unsigned xGlobal, const unsigned yGlobal,
                       const unsigned zGlobal) const;
    std::optional<Fields> valueAtUp(const double xGlobal, const double yGlobal,
                                    const double zGlobal) const;
    void removeChildren();
    void addParent(AMRNode<SolverType, ProblemType, Fields, padding>* parent);
    void addSibling(AMRNode<SolverType, ProblemType, Fields, padding>* sibling);
    void addChild(AMRNode<SolverType, ProblemType, Fields, padding>* child);

    PaddedGrid<FieldStruct, GHOST_CELLS>& grid;
    GridBoundary::Boundary gridBoundary;

    Solver<SolverType, Fields, ProblemType, padding> solver;

  private:
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> parents;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> siblings;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> children;

    const AMRParameters& configuration;

    Fields interpolateDown(const unsigned xLocal, const unsigned yLocal,
                           const unsigned zLocal) const;
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
void AMRNode<SolverType, ProblemType, Fields, padding>::setValue(const unsigned xGlobal,
                                                                 const unsigned yGlobal,
                                                                 const unsigned zGlobal,
                                                                 const Fields& fields) {
    const unsigned x = padding + xGlobal - this->gridBoundary[0].first;
    const unsigned y = padding + yGlobal - this->gridBoundary[1].first;
    const unsigned z = padding + zGlobal - this->gridBoundary[2].first;
    this->grid(x, y, z) = fields;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::injectParents() {
    GridBoundary::Boundary translated =
        GridBoundary::translateDown(this->gridBoundary, this->configuration.refinementFactor);
    for (AMRNode<SolverType, ProblemType, Fields, padding>* parent : this->parents) {
        std::optional<GridBoundary::Boundary> overlappingArea =
            GridBoundary::overlappingArea(translated, parent->gridBoundary);
        if (overlappingArea.has_value()) {
            GridBoundary::Boundary& overlap = overlappingArea.value();
            for (unsigned x = overlap[0].first; x <= overlap[0].second; x++) {
                for (unsigned y = overlap[1].first; y <= overlap[1].second; y++) {
                    for (unsigned z = overlap[2].first; z <= overlap[2].second; z++) {
                        parent->setValue(x, y, z, this->valueAtDown(x, y, z));
                    }
                }
            }
        }
    }
}

/* Careful, this function does not check whether the requested indices are valid. */
template <class SolverType, class ProblemType, class Fields, unsigned padding>
Fields AMRNode<SolverType, ProblemType, Fields, padding>::valueAtDown(
    const unsigned xGlobal, const unsigned yGlobal, const unsigned zGlobal) const {
    return this->interpolateDown(xGlobal * this->configuration.refinementFactor,
                                 yGlobal * this->configuration.refinementFactor,
                                 zGlobal * this->configuration.refinementFactor);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
Fields AMRNode<SolverType, ProblemType, Fields, padding>::interpolateDown(
    const unsigned xLocal, const unsigned yLocal, const unsigned zLocal) const {
    const unsigned xStart = std::max(xLocal, this->gridBoundary[0].first);
    const unsigned xEnd =
        std::min(xLocal + this->configuration.refinementFactor, this->gridBoundary[0].second + 1);
    const unsigned yStart = std::max(yLocal, this->gridBoundary[1].first);
    const unsigned yEnd =
        std::min(yLocal + this->configuration.refinementFactor, this->gridBoundary[1].second + 1);
    const unsigned zStart = std::max(zLocal, this->gridBoundary[2].first);
    const unsigned zEnd =
        std::min(zLocal + this->configuration.refinementFactor, this->gridBoundary[2].second + 1);
    Fields fields = { 0.0 };
    for (unsigned x = xStart; x < xEnd; x++) {
        for (unsigned y = yStart; y < yEnd; y++) {
            for (unsigned z = zStart; z < zEnd; z++) {
                for (unsigned field = 0; field < fields.size(); field++) {
                    fields[field] += this->grid(x, y, z)[field];
                }
            }
        }
    }
    const double normalisationFactor =
        1.0 / (static_cast<double>(xEnd - xStart) * static_cast<double>(yEnd - yStart) *
               static_cast<double>(zEnd - zStart));
    for (unsigned field = 0; field < fields.size(); field++) {
        fields[field] *= normalisationFactor;
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
