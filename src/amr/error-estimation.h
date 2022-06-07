#pragma once

#include <optional>

#include "../boundary/boundary.h"
#include "../grid/padded-grid.h"
#include "../solver/runge-kutta-solver.h"
#include "flags.h"

namespace ErrorEstimation {
template <class ProblemType, class Fields, unsigned padding>
std::optional<CellFlags::Flags>
getFlags(PaddedGrid<Fields, padding>& grid, const std::array<unsigned, 3> gridOffset,
         const unsigned timeDelta, const double errorThreshold, const Problem<ProblemType>& problem,
         const std::array<BoundaryType, Faces::FaceMax>& boundaryTypes);
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

template <class Fields, unsigned padding>
PaddedGrid<Fields, padding>
createUpper(const PaddedGrid<Fields, padding>& grid,
            const std::array<BoundaryType, Faces::FaceMax> boundaryTypes) {
    PaddedGrid<Fields, padding> upper(
        grid.defaultValue, (grid.xEnd() - grid.xStart()) / 2, (grid.yEnd() - grid.yStart()) / 2,
        (grid.zEnd() - grid.zStart()) / 2, grid.posLeft, grid.posRight, boundaryTypes);
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
    return upper;
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
ErrorEstimation::getFlags(PaddedGrid<Fields, padding>& grid,
                          const std::array<unsigned, 3> gridOffset, const unsigned timeDelta,
                          const double errorThreshold, const Problem<ProblemType>& problem,
                          const std::array<BoundaryType, Faces::FaceMax>& boundaryTypes) {
    PaddedGrid<Fields, padding> baseGrid(grid);
    PaddedGrid<Fields, padding> upperGrid = createUpper(baseGrid, boundaryTypes);
    RungeKuttaSolver<ProblemType, Fields, padding> solver1(baseGrid, problem, 1);
    solver1.timeDelta = timeDelta;
    solver1.singleStep();
    solver1.singleStep();
    RungeKuttaSolver<ProblemType, Fields, padding> solver2(upperGrid, problem, 1);
    solver2.timeDelta = timeDelta * 2.0;
    solver2.singleStep();

    CellFlags::FlagsX flagsX;
    unsigned firstX = 0;
    bool foundError = false;
    for (unsigned x = upperGrid.xStart(); x < upperGrid.xEnd(); x++) {
        CellFlags::FlagsY flagsY;
        unsigned firstY = 0;
        bool foundY = false;
        for (unsigned y = upperGrid.yStart(); y < upperGrid.yEnd(); y++) {
            CellFlags::FlagsZ flagsZ;
            unsigned start = 0;
            bool streak = false;
            for (unsigned z = upperGrid.zStart(); z < upperGrid.zEnd(); z++) {
                if (checkThreshold(upperGrid, baseGrid, x, y, z, errorThreshold)) {
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
                flagsZ.push_back(std::make_pair(2 * start - padding + gridOffset[2],
                                                2 * (grid.zEnd() - 1) - padding + gridOffset[2]));
            }
            if (foundY) {
                flagsY.push_back(flagsZ);
            }
        }
        if (foundError) {
            flagsX.push_back(std::make_pair(firstY, flagsY));
        }
    }
    if (foundError) {
        return std::make_optional(std::make_pair(firstX, flagsX));
    } else {
        return std::nullopt;
    }
}
