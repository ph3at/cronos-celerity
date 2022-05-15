#pragma once

#include "../data-types/faces.h"
#include "../data-types/phys-fields.h"
#include "../grid/padded-grid.h"

namespace Outflow {
template <unsigned padding>
void apply(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const unsigned face);

template <unsigned padding>
void applyX(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);

template <unsigned padding>
void applyY(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);

template <unsigned padding>
void applyZ(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft);
}; // namespace Outflow

template <unsigned padding>
void Outflow::applyX(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                     const bool isLeft) {
    const bool isParallelVelocity = field == FieldNames::VELOCITY_X;
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
                if (grid(xInner - direction, y, z)[FieldNames::VELOCITY_X] * direction > 0.0) {
                    grid(x, y, z)[field] = grid(x - direction, y, z)[field];
                } else if (isParallelVelocity) {
                    grid(x, y, z)[field] = -grid(2 * xInner - x - direction, y, z)[field];
                } else {
                    grid(x, y, z)[field] = grid(2 * xInner - x - direction, y, z)[field];
                }
            }
        }
    }
}

template <unsigned padding>
void Outflow::applyY(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                     const bool isLeft) {
    const bool isParallelVelocity = field == FieldNames::VELOCITY_Y;
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
                if (grid(x, yInner - direction, z)[FieldNames::VELOCITY_Y] * direction > 0.0) {
                    grid(x, y, z)[field] = grid(x, y - direction, z)[field];
                } else if (isParallelVelocity) {
                    grid(x, y, z)[field] = -grid(x, 2 * yInner - y - direction, z)[field];
                } else {
                    grid(x, y, z)[field] = grid(x, 2 * yInner - y - direction, z)[field];
                }
            }
        }
    }
}

template <unsigned padding>
void Outflow::applyZ(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                     const bool isLeft) {
    const bool isParallelVelocity = field == FieldNames::VELOCITY_Z;
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
                if (grid(x, y, zInner - direction)[FieldNames::VELOCITY_Z] * direction > 0.0) {
                    grid(x, y, z)[field] = grid(x, y, z - direction)[field];
                } else if (isParallelVelocity) {
                    grid(x, y, z)[field] = -grid(x, y, 2 * zInner - z - direction)[field];
                } else {
                    grid(x, y, z)[field] = grid(x, y, 2 * zInner - z - direction)[field];
                }
            }
        }
    }
}

template <unsigned padding>
void Outflow::apply(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
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

