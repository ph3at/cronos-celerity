#pragma once

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
inline Fields interpolate(const PaddedGrid<Fields, padding>& other, unsigned x, unsigned y,
                          unsigned z) {
    const unsigned xBase = 2 * x;
    const unsigned yBase = 2 * y;
    const unsigned zBase = 2 * z;
    Fields fields = { 0.0 };
    for (unsigned field = 0; field < fields.size(); field++) {
        fields[field] += other(xBase, yBase, zBase)[field];
        fields[field] += other(xBase, yBase, zBase + 1)[field];
        fields[field] += other(xBase, yBase + 1, zBase)[field];
        fields[field] += other(xBase, yBase + 1, zBase + 1)[field];
        fields[field] += other(xBase + 1, yBase, zBase)[field];
        fields[field] += other(xBase + 1, yBase, zBase + 1)[field];
        fields[field] += other(xBase + 1, yBase + 1, zBase)[field];
        fields[field] += other(xBase + 1, yBase + 1, zBase + 1)[field];
        fields[field] /= 8.0;
    }
    return fields;
}

template <class ProblemType> ProblemType createUpperProblem(const ProblemType& problem) {
    ProblemType upperProblem(problem);
    for (unsigned dir = 0; dir < Direction::DirMax; dir++) {
        upperProblem.numberCells[dir] = problem.numberCells[dir] / 2;
        upperProblem.cellSize[dir] = (upperProblem.posRight[dir] - upperProblem.posLeft[dir]) /
                                     upperProblem.numberCells[dir];
        upperProblem.inverseCellSize[dir] = 1.0 / upperProblem.cellSize[dir];
    }
    return upperProblem;
}

template <class Fields, unsigned padding>
void initialiseUpperGrid(PaddedGrid<Fields, padding>& upper,
                         const PaddedGrid<Fields, padding>& grid) {
    for (unsigned x = 1; x + 1 < upper.xDim(); x++) {
        for (unsigned y = 1; y + 1 < upper.yDim(); y++) {
            for (unsigned z = 1; z + 1 < upper.zDim(); z++) {
                upper(x, y, z) = interpolate(grid, x - 1, y - 1, z - 1);
            }
            upper(x, y, 0) = upper(x, y, 1);
            upper(x, y, upper.zDim() - 1) = upper(x, y, upper.zDim() - 2);
        }
        for (unsigned z = 0; z < upper.zDim(); z++) {
            upper(x, 0, z) = upper(x, 1, z);
            upper(x, upper.yDim() - 1, z) = upper(x, upper.yDim() - 2, z);
        }
    }
    for (unsigned y = 1; y + 1 < upper.yDim(); y++) {
        for (unsigned z = 1; z + 1 < upper.zDim(); z++) {
            upper(0, y, z) = upper(1, y, z);
            upper(upper.xDim() - 1, y, z) = upper(upper.xDim() - 2, y, z);
        }
    }
}

template <class Fields, unsigned padding>
inline bool checkThreshold(const PaddedGrid<Fields, padding>& upper,
                           const PaddedGrid<Fields, padding>& base, const unsigned x,
                           const unsigned y, const unsigned z, const double errorThreshold) {
    for (unsigned xx = 2 * x - padding; xx < 2 * x - padding + 2; xx++) {
        for (unsigned yy = 2 * y - padding; yy < 2 * y - padding + 2; yy++) {
            for (unsigned zz = 2 * z - padding; zz < 2 * z - padding + 2; zz++) {
                for (unsigned field = 0; field < Fields().size(); field++) {
                    const double diff = upper(x, y, z)[field] - base(xx, yy, zz)[field];
                    const double avg = (upper(x, y, z)[field] + base(xx, yy, zz)[field]) / 2.0;
                    double error;
                    if (avg == 0.0) {
                        error = std::abs(diff);
                    } else {
                        error = std::abs(diff / avg);
                    }
                    if (error > errorThreshold) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

template <class ProblemType, class Fields, unsigned padding>
std::optional<CellFlags::Flags>
ErrorEstimation::getFlags(const PaddedGrid<Fields, padding>& grid,
                          const std::array<unsigned, 3> gridOffset, const double timeDelta,
                          const double errorThreshold, const ProblemType& problem) {
    RungeKuttaSolver<ProblemType, Fields, padding> solver1(grid, problem, 1);
    solver1.timeDelta = timeDelta;
    solver1.singleStep();
    solver1.singleStep();

    ProblemType upperProblem = createUpperProblem(problem);
    RungeKuttaSolver<ProblemType, Fields, padding> solver2(upperProblem, 1);
    initialiseUpperGrid(solver2.grid, grid);
    solver2.timeDelta = timeDelta * 2.0;
    solver2.singleStep();

    CellFlags::FlagsX flagsX;
    unsigned firstX = 0;
    bool foundError = false;
    unsigned emptyYs = 0;
    for (unsigned x = solver2.grid.xStart(); x < solver2.grid.xEnd(); x++) {
        CellFlags::FlagsY flagsY;
        unsigned firstY = 0;
        bool foundY = false;
        unsigned emptyZs = 0;
        for (unsigned y = solver2.grid.yStart(); y < solver2.grid.yEnd(); y++) {
            CellFlags::FlagsZ flagsZ;
            unsigned start = 0;
            bool streak = false;
            for (unsigned z = solver2.grid.zStart(); z < solver2.grid.zEnd(); z++) {
                if (checkThreshold(solver2.grid, solver1.grid, x, y, z, errorThreshold)) {
                    if (!foundError) {
                        firstX = 2 * x - padding + gridOffset[0];
                        foundError = true;
                    }
                    if (!foundY) {
                        firstY = 2 * y - padding + gridOffset[1];
                        foundY = true;
                    }
                    if (!streak) {
                        start = z;
                        streak = true;
                    }
                } else if (streak) {
                    streak = false;
                    flagsZ.push_back(std::make_pair(2 * start - padding + gridOffset[2],
                                                    2 * z - padding + gridOffset[2]));
                }
            }
            if (streak) {
                flagsZ.push_back(
                    std::make_pair(2 * start - padding + gridOffset[2],
                                   2 * (solver2.grid.zEnd() - 1) - padding + gridOffset[2]));
            }
            if (foundY) {
                if (flagsZ.size() > 0) {
                    for (unsigned z = 0; z < emptyZs; z++) {
                        flagsY.push_back(CellFlags::FlagsZ({}));
                    }
                    emptyZs = 0;
                    flagsY.push_back(flagsZ);
                } else {
                    emptyZs++;
                }
            }
        }
        if (foundError) {
            if (foundY) {
                for (unsigned y = 0; y < emptyYs; y++) {
                    flagsX.push_back(std::make_pair(0, CellFlags::FlagsY({})));
                }
                emptyYs = 0;
                flagsX.push_back(std::make_pair(firstY, flagsY));
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
