#pragma once

#include <cmath>
#include <optional>

#include "../boundary/boundary.h"
#include "../grid/padded-grid.h"
#include "../solver/runge-kutta-solver.h"
#include "flags.h"

namespace ErrorEstimation {
template <class ProblemType, class Fields, unsigned padding>
std::optional<CellFlags::Flags>
getFlags(const PaddedGrid<Fields, padding>& grid, const std::array<unsigned, 3> gridOffset,
         const double timeDelta, const double errorThreshold, const ProblemType& problem);
};

template <class Fields, unsigned padding>
inline Fields average(const PaddedGrid<Fields, padding>& baseGrid, unsigned x, unsigned y,
                      unsigned z) {
    const size_t xBase = 2 * x - padding;
    const size_t yBase = 2 * y - padding;
    const size_t zBase = 2 * z - padding;
    Fields fields = { 0.0 };
    double lowerCells = 0.0;
    for (size_t xGrid = xBase; xGrid < std::min(xBase + 2, baseGrid.xEnd()); xGrid++) {
        for (size_t yGrid = yBase; yGrid < std::min(yBase + 2, baseGrid.yEnd()); yGrid++) {
            for (size_t zGrid = zBase; zGrid < std::min(zBase + 2, baseGrid.zEnd()); zGrid++) {
                for (size_t field = 0; field < fields.size(); field++) {
                    fields[field] += baseGrid(xGrid, yGrid, zGrid)[field];
                }
                lowerCells += 1.0;
            }
        }
    }
    for (unsigned field = 0; field < fields.size(); field++) {
        fields[field] /= lowerCells;
    }
    return fields;
}

template <class ProblemType> ProblemType createUpperProblem(const ProblemType& problem) {
    ProblemType upperProblem(problem);
    for (unsigned dir = 0; dir < Direction::DirMax; dir++) {
        upperProblem.numberCells[dir] = (problem.numberCells[dir] + 1) / 2;
        upperProblem.cellSize[dir] = problem.cellSize[dir] * 2.0;
        upperProblem.inverseCellSize[dir] = problem.inverseCellSize[dir] / 2.0;
        upperProblem.posRight[dir] =
            upperProblem.posLeft[dir] + upperProblem.numberCells[dir] * upperProblem.cellSize[dir];
    }
    return upperProblem;
}

template <class Fields, unsigned padding>
void initialiseUpperGrid(PaddedGrid<Fields, padding>& upper,
                         const PaddedGrid<Fields, padding>& grid) {
    for (unsigned x = upper.xStart(); x < upper.xEnd(); x++) {
        for (unsigned y = upper.yStart(); y < upper.yEnd(); y++) {
            for (unsigned z = upper.zStart(); z < upper.zEnd(); z++) {
                upper(x, y, z) = average(grid, x, y, z);
            }
        }
    }
}

template <class Fields, unsigned padding>
inline bool checkThreshold(const PaddedGrid<Fields, padding>& upper,
                           const PaddedGrid<Fields, padding>& lower, const unsigned x,
                           const unsigned y, const unsigned z, const double errorThreshold) {
    Fields lowerValues = average(lower, x, y, z);
    for (unsigned field = 0; field < Fields().size(); field++) {
        const double diff = std::abs(upper(x, y, z)[field] - lowerValues[field]);
        const double avg = 0.5 * (std::abs(upper(x, y, z)[field]) + std::abs(lowerValues[field]));
        double error;
        if (avg != 0.0) {
            error = diff / avg;
        } else if (lowerValues[field] != 0.0) {
            error = diff / lowerValues[field];
        } else {
            error = diff;
        }
        if (error > errorThreshold) {
            return true;
        }
    }
    return false;
}

template <class ProblemType, class Fields, unsigned padding>
std::optional<CellFlags::Flags>
ErrorEstimation::getFlags(const PaddedGrid<Fields, padding>& grid,
                          const std::array<unsigned, 3> gridOffset, const double timeDelta,
                          const double errorThreshold, const ProblemType& problem) {
    RungeKuttaSolver<ProblemType, Fields, padding> solverRegular(grid, problem, 1);
    Boundary::applyAll(solverRegular.grid, problem);
    solverRegular.timeDelta = timeDelta;
    solverRegular.singleStep();
    solverRegular.singleStep();

    ProblemType upperProblem = createUpperProblem(problem);
    RungeKuttaSolver<ProblemType, Fields, padding> solverCoarse(upperProblem, 1);
    initialiseUpperGrid(solverCoarse.grid, grid);
    Boundary::applyAll(solverCoarse.grid, upperProblem);
    solverCoarse.timeDelta = timeDelta * 2.0;
    solverCoarse.singleStep();

    const unsigned xMax = gridOffset[Direction::DirX] + grid.xEnd() - padding - 1;
    const unsigned yMax = gridOffset[Direction::DirY] + grid.yEnd() - padding - 1;
    const unsigned zMax = gridOffset[Direction::DirZ] + grid.zEnd() - padding - 1;

    CellFlags::FlagsX flagsX;
    unsigned firstX = 0;
    bool foundError = false;
    unsigned emptyYs = 0;
    for (unsigned x = solverCoarse.grid.xStart(); x < solverCoarse.grid.xEnd(); x++) {
        CellFlags::FlagsY flagsY;
        unsigned firstY = 0;
        bool foundY = false;
        unsigned emptyZs = 0;
        for (unsigned y = solverCoarse.grid.yStart(); y < solverCoarse.grid.yEnd(); y++) {
            CellFlags::FlagsZ flagsZ;
            unsigned start = 0;
            bool streak = false;
            for (unsigned z = solverCoarse.grid.zStart(); z < solverCoarse.grid.zEnd(); z++) {
                if (checkThreshold(solverCoarse.grid, solverRegular.grid, x, y, z,
                                   errorThreshold)) {
                    if (!foundError) {
                        firstX = 2 * (x - padding) + gridOffset[0];
                        foundError = true;
                    }
                    if (!foundY) {
                        firstY = 2 * (y - padding) + gridOffset[1];
                        foundY = true;
                    }
                    if (!streak) {
                        start = z;
                        streak = true;
                    }
                } else if (streak) {
                    streak = false;
                    flagsZ.push_back(std::make_pair(2 * (start - padding) + gridOffset[2],
                                                    2 * (z - padding) + gridOffset[2] - 1));
                }
            }
            if (streak) {
                flagsZ.push_back(std::make_pair(2 * (start - padding) + gridOffset[2], zMax));
            }
            if (foundY) {
                if (flagsZ.size() > 0) {
                    for (unsigned z = 0; z < emptyZs; z++) {
                        flagsY.emplace_back();
                        flagsY.emplace_back();
                    }
                    emptyZs = 0;
                    flagsY.push_back(flagsZ);
                    if (firstY + flagsY.size() <= yMax) {
                        flagsY.push_back(flagsZ);
                    }
                } else {
                    emptyZs++;
                }
            }
        }
        if (foundError) {
            if (foundY) {
                for (unsigned y = 0; y < emptyYs; y++) {
                    flagsX.emplace_back(0, CellFlags::FlagsY());
                    flagsX.emplace_back(0, CellFlags::FlagsY());
                }
                emptyYs = 0;
                flagsX.push_back(std::make_pair(firstY, flagsY));
                if (firstX + flagsX.size() <= xMax) {
                    flagsX.push_back(std::make_pair(firstY, flagsY));
                }
            } else {
                emptyYs++;
            }
        }
    }
    if (foundError) {
        return std::make_optional(std::make_pair(firstX, flagsX));
    } else {
        return std::nullopt;
    }
}
