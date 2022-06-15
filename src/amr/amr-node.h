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
    AMRNode(const unsigned id, const PaddedGrid<Fields, padding>& grid,
            const GridBoundary::Boundary gridBoundary,
            const Problem<ProblemType, Fields, padding>& problem,
            const AMRParameters& configuration, const double timeDelta, const double timeCurrent);
    AMRNode(const unsigned id, const Fields defaultValue, const std::array<unsigned, 3> dim,
            const std::array<double, 3> posLeft, const std::array<double, 3> posRight,
            const GridBoundary::Boundary gridBoundary,
            const Problem<ProblemType, Fields, padding>& problem,
            const AMRParameters& configuration, const double timeDelta, const double timeCurrent);

    void setValue(const unsigned xGlobal, const unsigned yGlobal, const unsigned zGlobal,
                  const Fields& fields);
    void injectParents();
    void updateBoundary();
    void updateBoundarySiblings();
    Fields valueAtDown(const unsigned xGlobal, const unsigned yGlobal,
                       const unsigned zGlobal) const;
    std::optional<Fields> valueAtUp(const double xGlobal, const double yGlobal,
                                    const double zGlobal) const;
    void removeChildren();
    void addParent(AMRNode<SolverType, ProblemType, Fields, padding>* parent);
    void addSibling(AMRNode<SolverType, ProblemType, Fields, padding>* sibling);
    void addChild(AMRNode<SolverType, ProblemType, Fields, padding>* child);

    const unsigned id;
    PaddedGrid<FieldStruct, GHOST_CELLS> grid;
    GridBoundary::Boundary gridBoundary;

    SolverType solver;

    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> parents;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> siblings;
    std::vector<AMRNode<SolverType, ProblemType, Fields, padding>*> children;

  private:
    const AMRParameters& configuration;

    Fields interpolateDown(const unsigned xLocal, const unsigned yLocal,
                           const unsigned zLocal) const;
    Fields interpolateUp(const double xLocal, const double yLocal, const double zLocal) const;

    void boundaryXLow(const bool includeParents);
    void boundaryXHigh(const bool includeParents);
    void boundaryYLow(const bool includeParents);
    void boundaryYHigh(const bool includeParents);
    void boundaryZLow(const bool includeParents);
    void boundaryZHigh(const bool includeParents);
};

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>::AMRNode(
    const unsigned id, const PaddedGrid<Fields, padding>& grid,
    const GridBoundary::Boundary gridBoundary, const Problem<ProblemType, Fields, padding>& problem,
    const AMRParameters& configuration, const double timeDelta, const double timeCurrent)
    : id(id), grid(grid), gridBoundary(gridBoundary), solver(SolverType(this->grid, problem)),
      configuration(configuration) {
    this->solver.timeDelta = timeDelta;
    this->solver.timeCurrent = timeCurrent;
    this->solver.timeStep = 0;
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
AMRNode<SolverType, ProblemType, Fields, padding>::AMRNode(
    const unsigned id, const Fields defaultValue, const std::array<unsigned, 3> dim,
    const std::array<double, 3> posLeft, const std::array<double, 3> posRight,
    const GridBoundary::Boundary gridBoundary, const Problem<ProblemType, Fields, padding>& problem,
    const AMRParameters& configuration, const double timeDelta, const double timeCurrent)
    : id(id), grid(defaultValue, dim[0], dim[1], dim[2], posLeft, posRight),
      gridBoundary(gridBoundary), solver(SolverType(this->grid, problem)),
      configuration(configuration) {
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

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::updateBoundary() {
    const bool includeParents = true;
    this->boundaryXLow(includeParents);
    this->boundaryXHigh(includeParents);
    this->boundaryYLow(includeParents);
    this->boundaryYHigh(includeParents);
    this->boundaryZLow(includeParents);
    this->boundaryZHigh(includeParents);
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::updateBoundarySiblings() {
    const bool includeParents = false;
    this->boundaryXLow(includeParents);
    this->boundaryXHigh(includeParents);
    this->boundaryYLow(includeParents);
    this->boundaryYHigh(includeParents);
    this->boundaryZLow(includeParents);
    this->boundaryZHigh(includeParents);
}

/* Careful, this function does not check whether the requested indices are valid. */
template <class SolverType, class ProblemType, class Fields, unsigned padding>
Fields AMRNode<SolverType, ProblemType, Fields, padding>::valueAtDown(
    const unsigned xGlobal, const unsigned yGlobal, const unsigned zGlobal) const {
    return this->interpolateDown(
        padding + xGlobal * this->configuration.refinementFactor - this->gridBoundary[0].first,
        padding + yGlobal * this->configuration.refinementFactor - this->gridBoundary[1].first,
        padding + zGlobal * this->configuration.refinementFactor - this->gridBoundary[2].first);
}

/* Right now, for every cell in the lower grid is entirely filled out in the upper grid. If that
 * changes, one should use ints instead. */
template <class SolverType, class ProblemType, class Fields, unsigned padding>
Fields AMRNode<SolverType, ProblemType, Fields, padding>::interpolateDown(
    const unsigned xLocal, const unsigned yLocal, const unsigned zLocal) const {
    const unsigned xStart = std::max(xLocal, static_cast<unsigned>(this->grid.xStart()));
    const unsigned xEnd = std::min(xLocal + this->configuration.refinementFactor,
                                   static_cast<unsigned>(this->grid.xEnd()));
    const unsigned yStart = std::max(yLocal, static_cast<unsigned>(this->grid.yStart()));
    const unsigned yEnd = std::min(yLocal + this->configuration.refinementFactor,
                                   static_cast<unsigned>(this->grid.yEnd()));
    const unsigned zStart = std::max(zLocal, static_cast<unsigned>(this->grid.zStart()));
    const unsigned zEnd = std::min(zLocal + this->configuration.refinementFactor,
                                   static_cast<unsigned>(this->grid.zEnd()));
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
    const double xLocal =
        xGlobal - static_cast<double>(this->gridBoundary[0].first) + static_cast<double>(padding);
    const double yLocal =
        yGlobal - static_cast<double>(this->gridBoundary[1].first) + static_cast<double>(padding);
    const double zLocal =
        zGlobal - static_cast<double>(this->gridBoundary[2].first) + static_cast<double>(padding);
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

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::boundaryXLow(const bool includeParents) {
    for (int x = static_cast<int>(this->gridBoundary[0].first) - padding;
         x < static_cast<int>(this->gridBoundary[0].first); x++) {
        const unsigned xLocal = static_cast<unsigned>(x + padding) - this->gridBoundary[0].first;
        for (int y = static_cast<int>(this->gridBoundary[1].first) - padding;
             y <= static_cast<int>(this->gridBoundary[1].second + padding); y++) {
            const unsigned yLocal =
                static_cast<unsigned>(y + padding) - this->gridBoundary[1].first;
            for (int z = static_cast<int>(this->gridBoundary[2].first) - padding;
                 z <= static_cast<int>(this->gridBoundary[2].second + padding); z++) {
                unsigned zLocal = static_cast<unsigned>(z + padding) - this->gridBoundary[2].first;
                bool found = false;
                for (AMRNode<SolverType, ProblemType, Fields, padding>* sibling : this->siblings) {
                    if (x >= static_cast<int>(sibling->gridBoundary[0].first) &&
                        x <= static_cast<int>(sibling->gridBoundary[0].second) &&
                        y >= static_cast<int>(sibling->gridBoundary[1].first) &&
                        y <= static_cast<int>(sibling->gridBoundary[1].second) &&
                        z >= static_cast<int>(sibling->gridBoundary[2].first) &&
                        z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                        const unsigned xSibling =
                            static_cast<unsigned>(x + padding) - sibling->gridBoundary[0].first;
                        const unsigned ySibling =
                            static_cast<unsigned>(y + padding) - sibling->gridBoundary[1].first;
                        unsigned zSibling =
                            static_cast<unsigned>(z + padding) - sibling->gridBoundary[2].first;
                        while (z <= static_cast<int>(this->gridBoundary[2].second + padding) &&
                               z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                            this->grid(xLocal, yLocal, zLocal) =
                                sibling->grid(xSibling, ySibling, zSibling);
                            z++;
                            zLocal++;
                            zSibling++;
                        }
                        z--;
                        zLocal--;
                        found = true;
                        break;
                    }
                }
                if (includeParents && !found) {
                    const double xGlobal =
                        static_cast<double>(x) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double yGlobal =
                        static_cast<double>(y) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double zGlobal =
                        static_cast<double>(z) /
                        static_cast<double>(this->configuration.refinementFactor);
                    for (AMRNode<SolverType, ProblemType, Fields, padding>* parent :
                         this->parents) {
                        std::optional<Fields> fields = parent->valueAtUp(xGlobal, yGlobal, zGlobal);
                        if (fields.has_value()) {
                            this->grid(xLocal, yLocal, zLocal) = fields.value();
                            break;
                        }
                    }
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::boundaryXHigh(const bool includeParents) {
    for (int x = static_cast<int>(this->gridBoundary[0].second) + 1;
         x <= static_cast<int>(this->gridBoundary[0].second + padding); x++) {
        const unsigned xLocal = static_cast<unsigned>(x + padding) - this->gridBoundary[0].first;
        for (int y = static_cast<int>(this->gridBoundary[1].first) - padding;
             y <= static_cast<int>(this->gridBoundary[1].second + padding); y++) {
            const unsigned yLocal =
                static_cast<unsigned>(y + padding) - this->gridBoundary[1].first;
            for (int z = static_cast<int>(this->gridBoundary[2].first) - padding;
                 z <= static_cast<int>(this->gridBoundary[2].second + padding); z++) {
                unsigned zLocal = static_cast<unsigned>(z + padding) - this->gridBoundary[2].first;
                bool found = false;
                for (AMRNode<SolverType, ProblemType, Fields, padding>* sibling : this->siblings) {
                    if (x >= static_cast<int>(sibling->gridBoundary[0].first) &&
                        x <= static_cast<int>(sibling->gridBoundary[0].second) &&
                        y >= static_cast<int>(sibling->gridBoundary[1].first) &&
                        y <= static_cast<int>(sibling->gridBoundary[1].second) &&
                        z >= static_cast<int>(sibling->gridBoundary[2].first) &&
                        z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                        const unsigned xSibling =
                            static_cast<unsigned>(x + padding) - sibling->gridBoundary[0].first;
                        const unsigned ySibling =
                            static_cast<unsigned>(y + padding) - sibling->gridBoundary[1].first;
                        unsigned zSibling =
                            static_cast<unsigned>(z + padding) - sibling->gridBoundary[2].first;
                        while (z <= static_cast<int>(this->gridBoundary[2].second + padding) &&
                               z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                            this->grid(xLocal, yLocal, zLocal) =
                                sibling->grid(xSibling, ySibling, zSibling);
                            z++;
                            zLocal++;
                            zSibling++;
                        }
                        z--;
                        zLocal--;
                        found = true;
                        break;
                    }
                }
                if (includeParents && !found) {
                    const double xGlobal =
                        static_cast<double>(x) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double yGlobal =
                        static_cast<double>(y) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double zGlobal =
                        static_cast<double>(z) /
                        static_cast<double>(this->configuration.refinementFactor);
                    for (AMRNode<SolverType, ProblemType, Fields, padding>* parent :
                         this->parents) {
                        std::optional<Fields> fields = parent->valueAtUp(xGlobal, yGlobal, zGlobal);
                        if (fields.has_value()) {
                            this->grid(xLocal, yLocal, zLocal) = fields.value();
                            break;
                        }
                    }
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::boundaryYLow(const bool includeParents) {
    for (int y = static_cast<int>(this->gridBoundary[1].first) - padding;
         y < static_cast<int>(this->gridBoundary[1].first); y++) {
        const unsigned yLocal = static_cast<unsigned>(y + padding) - this->gridBoundary[1].first;
        for (int x = static_cast<int>(this->gridBoundary[0].first);
             x <= static_cast<int>(this->gridBoundary[0].second); x++) {
            const unsigned xLocal =
                static_cast<unsigned>(x + padding) - this->gridBoundary[0].first;
            for (int z = static_cast<int>(this->gridBoundary[2].first) - padding;
                 z <= static_cast<int>(this->gridBoundary[2].second + padding); z++) {
                unsigned zLocal = static_cast<unsigned>(z + padding) - this->gridBoundary[2].first;
                bool found = false;
                for (AMRNode<SolverType, ProblemType, Fields, padding>* sibling : this->siblings) {
                    if (x >= static_cast<int>(sibling->gridBoundary[0].first) &&
                        x <= static_cast<int>(sibling->gridBoundary[0].second) &&
                        y >= static_cast<int>(sibling->gridBoundary[1].first) &&
                        y <= static_cast<int>(sibling->gridBoundary[1].second) &&
                        z >= static_cast<int>(sibling->gridBoundary[2].first) &&
                        z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                        const unsigned xSibling =
                            static_cast<unsigned>(x + padding) - sibling->gridBoundary[0].first;
                        const unsigned ySibling =
                            static_cast<unsigned>(y + padding) - sibling->gridBoundary[1].first;
                        unsigned zSibling =
                            static_cast<unsigned>(z + padding) - sibling->gridBoundary[2].first;
                        while (z <= static_cast<int>(this->gridBoundary[2].second + padding) &&
                               z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                            this->grid(xLocal, yLocal, zLocal) =
                                sibling->grid(xSibling, ySibling, zSibling);
                            z++;
                            zLocal++;
                            zSibling++;
                        }
                        z--;
                        zLocal--;
                        found = true;
                        break;
                    }
                }
                if (includeParents && !found) {
                    const double xGlobal =
                        static_cast<double>(x) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double yGlobal =
                        static_cast<double>(y) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double zGlobal =
                        static_cast<double>(z) /
                        static_cast<double>(this->configuration.refinementFactor);
                    for (AMRNode<SolverType, ProblemType, Fields, padding>* parent :
                         this->parents) {
                        std::optional<Fields> fields = parent->valueAtUp(xGlobal, yGlobal, zGlobal);
                        if (fields.has_value()) {
                            this->grid(xLocal, yLocal, zLocal) = fields.value();
                            break;
                        }
                    }
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::boundaryYHigh(const bool includeParents) {
    for (int y = static_cast<int>(this->gridBoundary[1].second) + 1;
         y <= static_cast<int>(this->gridBoundary[1].second + padding); y++) {
        const unsigned yLocal = static_cast<unsigned>(y + padding) - this->gridBoundary[1].first;
        for (int x = static_cast<int>(this->gridBoundary[0].first);
             x <= static_cast<int>(this->gridBoundary[0].second); x++) {
            const unsigned xLocal =
                static_cast<unsigned>(x + padding) - this->gridBoundary[0].first;
            for (int z = static_cast<int>(this->gridBoundary[2].first) - padding;
                 z <= static_cast<int>(this->gridBoundary[2].second + padding); z++) {
                unsigned zLocal = static_cast<unsigned>(z + padding) - this->gridBoundary[2].first;
                bool found = false;
                for (AMRNode<SolverType, ProblemType, Fields, padding>* sibling : this->siblings) {
                    if (x >= static_cast<int>(sibling->gridBoundary[0].first) &&
                        x <= static_cast<int>(sibling->gridBoundary[0].second) &&
                        y >= static_cast<int>(sibling->gridBoundary[1].first) &&
                        y <= static_cast<int>(sibling->gridBoundary[1].second) &&
                        z >= static_cast<int>(sibling->gridBoundary[2].first) &&
                        z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                        const unsigned xSibling =
                            static_cast<unsigned>(x + padding) - sibling->gridBoundary[0].first;
                        const unsigned ySibling =
                            static_cast<unsigned>(y + padding) - sibling->gridBoundary[1].first;
                        unsigned zSibling =
                            static_cast<unsigned>(z + padding) - sibling->gridBoundary[2].first;
                        while (z <= static_cast<int>(this->gridBoundary[2].second + padding) &&
                               z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                            this->grid(xLocal, yLocal, zLocal) =
                                sibling->grid(xSibling, ySibling, zSibling);
                            z++;
                            zLocal++;
                            zSibling++;
                        }
                        z--;
                        zLocal--;
                        found = true;
                        break;
                    }
                }
                if (includeParents && !found) {
                    const double xGlobal =
                        static_cast<double>(x) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double yGlobal =
                        static_cast<double>(y) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double zGlobal =
                        static_cast<double>(z) /
                        static_cast<double>(this->configuration.refinementFactor);
                    for (AMRNode<SolverType, ProblemType, Fields, padding>* parent :
                         this->parents) {
                        std::optional<Fields> fields = parent->valueAtUp(xGlobal, yGlobal, zGlobal);
                        if (fields.has_value()) {
                            this->grid(xLocal, yLocal, zLocal) = fields.value();
                            break;
                        }
                    }
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::boundaryZLow(const bool includeParents) {
    for (int z = static_cast<int>(this->gridBoundary[2].first) - padding;
         z < static_cast<int>(this->gridBoundary[2].first); z++) {
        const unsigned zLocal = static_cast<unsigned>(z + padding) - this->gridBoundary[2].first;
        for (int x = static_cast<int>(this->gridBoundary[0].first);
             x <= static_cast<int>(this->gridBoundary[0].second); x++) {
            const unsigned xLocal =
                static_cast<unsigned>(x + padding) - this->gridBoundary[0].first;
            for (int y = static_cast<int>(this->gridBoundary[1].first);
                 y <= static_cast<int>(this->gridBoundary[1].second); y++) {
                unsigned yLocal = static_cast<unsigned>(y + padding) - this->gridBoundary[1].first;
                bool found = false;
                for (AMRNode<SolverType, ProblemType, Fields, padding>* sibling : this->siblings) {
                    if (x >= static_cast<int>(sibling->gridBoundary[0].first) &&
                        x <= static_cast<int>(sibling->gridBoundary[0].second) &&
                        y >= static_cast<int>(sibling->gridBoundary[1].first) &&
                        y <= static_cast<int>(sibling->gridBoundary[1].second) &&
                        z >= static_cast<int>(sibling->gridBoundary[2].first) &&
                        z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                        const unsigned xSibling =
                            static_cast<unsigned>(x + padding) - sibling->gridBoundary[0].first;
                        unsigned ySibling =
                            static_cast<unsigned>(y + padding) - sibling->gridBoundary[1].first;
                        const unsigned zSibling =
                            static_cast<unsigned>(z + padding) - sibling->gridBoundary[2].first;
                        while (y <= static_cast<int>(this->gridBoundary[1].second) &&
                               y <= static_cast<int>(sibling->gridBoundary[1].second)) {
                            this->grid(xLocal, yLocal, zLocal) =
                                sibling->grid(xSibling, ySibling, zSibling);
                            y++;
                            yLocal++;
                            ySibling++;
                        }
                        y--;
                        yLocal--;
                        found = true;
                        break;
                    }
                }
                if (includeParents && !found) {
                    const double xGlobal =
                        static_cast<double>(x) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double yGlobal =
                        static_cast<double>(y) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double zGlobal =
                        static_cast<double>(z) /
                        static_cast<double>(this->configuration.refinementFactor);
                    for (AMRNode<SolverType, ProblemType, Fields, padding>* parent :
                         this->parents) {
                        std::optional<Fields> fields = parent->valueAtUp(xGlobal, yGlobal, zGlobal);
                        if (fields.has_value()) {
                            this->grid(xLocal, yLocal, zLocal) = fields.value();
                            break;
                        }
                    }
                }
            }
        }
    }
}

template <class SolverType, class ProblemType, class Fields, unsigned padding>
void AMRNode<SolverType, ProblemType, Fields, padding>::boundaryZHigh(const bool includeParents) {
    for (int z = static_cast<int>(this->gridBoundary[2].second) + 1;
         z <= static_cast<int>(this->gridBoundary[2].second + padding); z++) {
        const unsigned zLocal = static_cast<unsigned>(z + padding) - this->gridBoundary[2].first;
        for (int x = static_cast<int>(this->gridBoundary[0].first);
             x <= static_cast<int>(this->gridBoundary[0].second); x++) {
            const unsigned xLocal =
                static_cast<unsigned>(x + padding) - this->gridBoundary[0].first;
            for (int y = static_cast<int>(this->gridBoundary[1].first);
                 y <= static_cast<int>(this->gridBoundary[1].second); y++) {
                unsigned yLocal = static_cast<unsigned>(y + padding) - this->gridBoundary[1].first;
                bool found = false;
                for (AMRNode<SolverType, ProblemType, Fields, padding>* sibling : this->siblings) {
                    if (x >= static_cast<int>(sibling->gridBoundary[0].first) &&
                        x <= static_cast<int>(sibling->gridBoundary[0].second) &&
                        y >= static_cast<int>(sibling->gridBoundary[1].first) &&
                        y <= static_cast<int>(sibling->gridBoundary[1].second) &&
                        z >= static_cast<int>(sibling->gridBoundary[2].first) &&
                        z <= static_cast<int>(sibling->gridBoundary[2].second)) {
                        const unsigned xSibling =
                            static_cast<unsigned>(x + padding) - sibling->gridBoundary[0].first;
                        unsigned ySibling =
                            static_cast<unsigned>(y + padding) - sibling->gridBoundary[1].first;
                        const unsigned zSibling =
                            static_cast<unsigned>(z + padding) - sibling->gridBoundary[2].first;
                        while (y <= static_cast<int>(this->gridBoundary[1].second) &&
                               y <= static_cast<int>(sibling->gridBoundary[1].second)) {
                            this->grid(xLocal, yLocal, zLocal) =
                                sibling->grid(xSibling, ySibling, zSibling);
                            y++;
                            yLocal++;
                            ySibling++;
                        }
                        y--;
                        yLocal--;
                        found = true;
                        break;
                    }
                }
                if (includeParents && !found) {
                    const double xGlobal =
                        static_cast<double>(x) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double yGlobal =
                        static_cast<double>(y) /
                        static_cast<double>(this->configuration.refinementFactor);
                    const double zGlobal =
                        static_cast<double>(z) /
                        static_cast<double>(this->configuration.refinementFactor);
                    for (AMRNode<SolverType, ProblemType, Fields, padding>* parent :
                         this->parents) {
                        std::optional<Fields> fields = parent->valueAtUp(xGlobal, yGlobal, zGlobal);
                        if (fields.has_value()) {
                            this->grid(xLocal, yLocal, zLocal) = fields.value();
                            break;
                        }
                    }
                }
            }
        }
    }
}
