#pragma once

#include <vector>

#include "../configuration/constants.h"
#include "../data-types/faces.h"
#include "../data-types/fields.h"
#include "../grid/padded-grid.h"
#include "../grid/utils.h"

namespace Reconstruction {

PerFaceValues reconstruct(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x, const unsigned y,
                          const unsigned z);

} // namespace Reconstruction

namespace ReconstructionSycl {

static inline double limitMinmod(const double deltaLeft, const double deltaRight) {
    if ((deltaLeft < 0 && deltaRight < 0) || (deltaLeft > 0 && deltaRight > 0)) {
        if (deltaRight > 0) {
            return std::min(deltaRight, deltaLeft);
        } else {
            return std::max(deltaRight, deltaLeft);
        }
    } else {
        return 0.0;
    }
}

static inline double derivX(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims, const unsigned x,
                            const unsigned y, const unsigned z, const unsigned field) {
    using grid::utils::idx3d;

    double deltaLeft = grid[idx3d(x, y, z, dims)][field] - grid[idx3d(x - 1, y, z, dims)][field];
    double deltaRight = grid[idx3d(x + 1, y, z, dims)][field] - grid[idx3d(x, y, z, dims)][field];

    return limitMinmod(deltaLeft, deltaRight);
}

static inline double derivY(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims, const unsigned x,
                            const unsigned y, const unsigned z, const unsigned field) {
    using grid::utils::idx3d;

    double deltaLeft = grid[idx3d(x, y, z, dims)][field] - grid[idx3d(x, y - 1, z, dims)][field];
    double deltaRight = grid[idx3d(x, y + 1, z, dims)][field] - grid[idx3d(x, y, z, dims)][field];

    return limitMinmod(deltaLeft, deltaRight);
}

static inline double derivZ(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims, const unsigned x,
                            const unsigned y, const unsigned z, const unsigned field) {
    using grid::utils::idx3d;

    double deltaLeft = grid[idx3d(x, y, z, dims)][field] - grid[idx3d(x, y, z - 1, dims)][field];
    double deltaRight = grid[idx3d(x, y, z + 1, dims)][field] - grid[idx3d(x, y, z, dims)][field];

    return limitMinmod(deltaLeft, deltaRight);
}

static inline PerFaceValues reconstruct(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
                                        const unsigned x, const unsigned y, const unsigned z) {
    using grid::utils::idx3d;

    PerFaceValues reconst = {};
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        reconst[Faces::FaceWest][field] = grid[idx3d(x, y, z, dims)][field] - 0.5 * derivX(grid, dims, x, y, z, field);
        reconst[Faces::FaceEast][field] =
            grid[idx3d(x - 1, y, z, dims)][field] + 0.5 * derivX(grid, dims, x - 1, y, z, field);

        reconst[Faces::FaceSouth][field] = grid[idx3d(x, y, z, dims)][field] - 0.5 * derivY(grid, dims, x, y, z, field);
        reconst[Faces::FaceNorth][field] =
            grid[idx3d(x, y - 1, z, dims)][field] + 0.5 * derivY(grid, dims, x, y - 1, z, field);

        reconst[Faces::FaceBottom][field] =
            grid[idx3d(x, y, z, dims)][field] - 0.5 * derivZ(grid, dims, x, y, z, field);
        reconst[Faces::FaceTop][field] =
            grid[idx3d(x, y, z - 1, dims)][field] + 0.5 * derivZ(grid, dims, x, y, z - 1, field);
    }

    return reconst;
}

} // namespace ReconstructionSycl
