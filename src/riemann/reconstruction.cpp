#include "reconstruction.h"

namespace Reconstruction {

double limitMinmod(const double deltaLeft, const double deltaRight) {
    if (deltaLeft * deltaRight > 0) {
        double deltaCenter = 0.5 * (deltaLeft + deltaRight);
        if (deltaRight > 0) {
            return std::min(deltaRight, std::min(deltaCenter, deltaLeft));
        } else {
            return std::max(deltaRight, std::max(deltaCenter, deltaLeft));
        }
    } else {
        return 0.0;
    }
}

double derivX(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x, const unsigned y, const unsigned z,
              const unsigned field) {
    double deltaLeft = grid(x, y, z)[field] - grid(x - 1, y, z)[field];
    double deltaRight = grid(x + 1, y, z)[field] - grid(x, y, z)[field];

    return limitMinmod(deltaLeft, deltaRight);
}

double derivY(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x, const unsigned y, const unsigned z,
              const unsigned field) {
    double deltaLeft = grid(x, y, z)[field] - grid(x, y - 1, z)[field];
    double deltaRight = grid(x, y + 1, z)[field] - grid(x, y, z)[field];

    return limitMinmod(deltaLeft, deltaRight);
}

double derivZ(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x, const unsigned y, const unsigned z,
              const unsigned field) {
    double deltaLeft = grid(x, y, z)[field] - grid(x, y, z - 1)[field];
    double deltaRight = grid(x, y, z + 1)[field] - grid(x, y, z)[field];

    return limitMinmod(deltaLeft, deltaRight);
}

PerFaceValues reconstruct(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const unsigned x, const unsigned y,
                          const unsigned z) {
    PerFaceValues reconst = {};
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        reconst[Faces::FaceWest][field] = grid(x, y, z)[field] - 0.5 * derivX(grid, x, y, z, field);
        reconst[Faces::FaceEast][field] = grid(x - 1, y, z)[field] + 0.5 * derivX(grid, x - 1, y, z, field);

        reconst[Faces::FaceSouth][field] = grid(x, y, z)[field] - 0.5 * derivY(grid, x, y, z, field);
        reconst[Faces::FaceNorth][field] = grid(x, y - 1, z)[field] + 0.5 * derivY(grid, x, y - 1, z, field);

        reconst[Faces::FaceBottom][field] = grid(x, y, z)[field] - 0.5 * derivZ(grid, x, y, z, field);
        reconst[Faces::FaceTop][field] = grid(x, y, z - 1)[field] + 0.5 * derivZ(grid, x, y, z - 1, field);
    }

    return reconst;
}

}; // namespace Reconstruction

namespace ReconstructionSycl {

double limitMinmod(const double deltaLeft, const double deltaRight) {
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

double derivX(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims, const unsigned x,
              const unsigned y, const unsigned z, const unsigned field) {
    using grid::utils::idx3d;

    double deltaLeft = grid[idx3d(x, y, z, dims)][field] - grid[idx3d(x - 1, y, z, dims)][field];
    double deltaRight = grid[idx3d(x + 1, y, z, dims)][field] - grid[idx3d(x, y, z, dims)][field];

    return limitMinmod(deltaLeft, deltaRight);
}

double derivY(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims, const unsigned x,
              const unsigned y, const unsigned z, const unsigned field) {
    using grid::utils::idx3d;

    double deltaLeft = grid[idx3d(x, y, z, dims)][field] - grid[idx3d(x, y - 1, z, dims)][field];
    double deltaRight = grid[idx3d(x, y + 1, z, dims)][field] - grid[idx3d(x, y, z, dims)][field];

    return limitMinmod(deltaLeft, deltaRight);
}

double derivZ(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims, const unsigned x,
              const unsigned y, const unsigned z, const unsigned field) {
    using grid::utils::idx3d;

    double deltaLeft = grid[idx3d(x, y, z, dims)][field] - grid[idx3d(x, y, z - 1, dims)][field];
    double deltaRight = grid[idx3d(x, y, z + 1, dims)][field] - grid[idx3d(x, y, z, dims)][field];

    return limitMinmod(deltaLeft, deltaRight);
}

PerFaceValues reconstruct(const std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims, const unsigned x,
                          const unsigned y, const unsigned z) {
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
