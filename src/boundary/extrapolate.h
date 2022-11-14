#pragma once

#include "../data-types/faces.h"
#include "../data-types/phys-fields.h"
#include "../grid/padded-grid.h"

namespace Extrapolate {
template <unsigned padding>
void apply(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const unsigned face);

template <unsigned padding>
void applyX(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);
template <unsigned padding>
void applyY(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);
template <unsigned padding>
void applyZ(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);
}; // namespace Extrapolate

template <unsigned padding>
void Extrapolate::applyX(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                         const bool isLeft) {
    int xInner, xOuter;
    int direction;
    if (isLeft) {
        xInner = grid.xStart() - 1;
        xOuter = -1;
        direction = -1;
    } else {
        xInner = grid.xEnd();
        xOuter = grid.xDim();
        direction = 1;
    }
    for (int x = xInner; x != xOuter; x += direction) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                grid(x, y, z)[field] = grid(x - direction, y, z)[field];
            }
        }
    }
}

template <unsigned padding>
void Extrapolate::applyY(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                         const bool isLeft) {
    int yInner, yOuter;
    int direction;
    if (isLeft) {
        yInner = grid.yStart() - 1;
        yOuter = -1;
        direction = -1;
    } else {
        yInner = grid.yEnd();
        yOuter = grid.yDim();
        direction = 1;
    }
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (int y = yInner; y != yOuter; y += direction) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                grid(x, y, z)[field] = grid(x, y - direction, z)[field];
            }
        }
    }
}

template <unsigned padding>
void Extrapolate::applyZ(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                         const bool isLeft) {
    int zInner, zOuter;
    int direction;
    if (isLeft) {
        zInner = grid.zStart() - 1;
        zOuter = -1;
        direction = -1;
    } else {
        zInner = grid.zEnd();
        zOuter = grid.zDim();
        direction = 1;
    }
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (int z = zInner; z != zOuter; z += direction) {
                grid(x, y, z)[field] = grid(x, y, z - direction)[field];
            }
        }
    }
}

template <unsigned padding>
void Extrapolate::apply(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                        const unsigned face) {
    if (face == Faces::FaceWest) {
        applyX(grid, field, true);
    } else if (face == Faces::FaceEast) {
        applyX(grid, field, false);
    } else if (face == Faces::FaceSouth) {
        applyY(grid, field, true);
    } else if (face == Faces::FaceNorth) {
        applyY(grid, field, false);
    } else if (face == Faces::FaceBottom) {
        applyZ(grid, field, true);
    } else {
        applyZ(grid, field, false);
    }
}

namespace ExtrapolateSycl {

template <unsigned padding>
void applyX(std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
            const unsigned field, const bool isLeft) {
    using grid::utils::idx3d;

    int xInner, xOuter;
    int direction;
    if (isLeft) {
        xInner = padding - 1;
        xOuter = -1;
        direction = -1;
    } else {
        xInner = dims[0] - padding;
        xOuter = dims[0];
        direction = 1;
    }
    for (int x = xInner; x != xOuter; x += direction) {
        for (unsigned y = 0; y < dims[1]; y++) {
            for (unsigned z = 0; z < dims[2]; z++) {
                grid[idx3d(x, y, z, dims)][field] = grid[idx3d(x - direction, y, z, dims)][field];
            }
        }
    }
}

template <unsigned padding>
void applyY(std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
            const unsigned field, const bool isLeft) {
    using grid::utils::idx3d;

    int yInner, yOuter;
    int direction;
    if (isLeft) {
        yInner = padding - 1;
        yOuter = -1;
        direction = -1;
    } else {
        yInner = dims[1] - padding;
        yOuter = dims[1];
        direction = 1;
    }
    for (unsigned x = 0; x < dims[0]; x++) {
        for (int y = yInner; y != yOuter; y += direction) {
            for (unsigned z = 0; z < dims[2]; z++) {
                grid[idx3d(x, y, z, dims)][field] = grid[idx3d(x, y - direction, z, dims)][field];
            }
        }
    }
}

template <unsigned padding>
void applyZ(std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
            const unsigned field, const bool isLeft) {
    using grid::utils::idx3d;

    int zInner, zOuter;
    int direction;
    if (isLeft) {
        zInner = padding - 1;
        zOuter = -1;
        direction = -1;
    } else {
        zInner = dims[2] - padding;
        zOuter = dims[2];
        direction = 1;
    }
    for (unsigned x = 0; x < dims[0]; x++) {
        for (unsigned y = 0; y < dims[1]; y++) {
            for (int z = zInner; z != zOuter; z += direction) {
                grid[idx3d(x, y, z, dims)][field] = grid[idx3d(x, y, z - direction, dims)][field];
            }
        }
    }
}

template <unsigned padding>
void apply(std::vector<FieldStruct>& grid, const grid::utils::dimensions& dims,
           const unsigned field, const unsigned face) {
    if (face == Faces::FaceWest) {
        applyX<padding>(grid, dims, field, true);
    } else if (face == Faces::FaceEast) {
        applyX<padding>(grid, dims, field, false);
    } else if (face == Faces::FaceSouth) {
        applyY<padding>(grid, dims, field, true);
    } else if (face == Faces::FaceNorth) {
        applyY<padding>(grid, dims, field, false);
    } else if (face == Faces::FaceBottom) {
        applyZ<padding>(grid, dims, field, true);
    } else {
        applyZ<padding>(grid, dims, field, false);
    }
}

} // namespace ExtrapolateSycl
