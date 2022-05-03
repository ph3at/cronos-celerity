#include "outflow.h"

#include "../misc/faces.h"

namespace Outflow {
template <unsigned padding>
void applyX(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft) {
    const bool isParallelVelocity = field == FieldNames::VELOCITY_X;
    unsigned xInner, xOuter;
    int direction;
    if (isLeft) {
        xInner = grid.xStart() - 1;
        xOuter = 0;
        direction = -1;
    } else {
        xInner = grid.xEnd();
        xOuter = grid.xDim() - 1;
        direction = 1;
    }
    for (unsigned x = xInner; x != xOuter; x += direction) {
        for (unsigned y = grid.yStart(); y < grid.yEnd(); y++) {
            for (unsigned z = grid.zStart(); z < grid.zEnd(); z++) {
                if (grid(xInner - direction, y, z)[FieldNames::VELOCITY_X] * direction > 0.0) {
                    grid(x, y, z)[field] = grid(x - direction, y, z)[field];
                } else if (isParallelVelocity) {
                    grid(x, y, z)[field] = -grid(2 * xInner - x - direction, y, z);
                } else {
                    grid(x, y, z)[field] = grid(2 * xInner - x - direction, y, z);
                }
            }
        }
    }
}

template <unsigned padding>
void applyY(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft) {
    const bool isParallelVelocity = field == FieldNames::VELOCITY_Y;
    unsigned yInner, yOuter;
    int direction;
    if (isLeft) {
        yInner = grid.yStart() - 1;
        yOuter = 0;
        direction = -1;
    } else {
        yInner = grid.yEnd();
        yOuter = grid.yDim() - 1;
        direction = 1;
    }
    for (unsigned x = grid.xStart(); x < grid.xEnd(); x++) {
        for (unsigned y = yInner; y != yOuter; y += direction) {
            for (unsigned z = grid.zStart(); z < grid.zEnd(); z++) {
                if (grid(x, yInner - direction, z)[FieldNames::VELOCITY_Y] * direction > 0.0) {
                    grid(x, y, z)[field] = grid(x, y - direction, z)[field];
                } else if (isParallelVelocity) {
                    grid(x, y, z)[field] = -grid(x, 2 * yInner - y - direction, z);
                } else {
                    grid(x, y, z)[field] = grid(x, 2 * yInner - y - direction, z);
                }
            }
        }
    }
}

template <unsigned padding>
void applyZ(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const bool isLeft) {
    const bool isParallelVelocity = field == FieldNames::VELOCITY_Z;
    unsigned zInner, zOuter;
    int direction;
    if (isLeft) {
        zInner = grid.zStart() - 1;
        zOuter = 0;
        direction = -1;
    } else {
        zInner = grid.zEnd();
        zOuter = grid.zDim() - 1;
        direction = 1;
    }
    for (unsigned x = grid.xStart(); x < grid.xEnd(); x++) {
        for (unsigned y = grid.yStart(); y < grid.yEnd(); y++) {
            for (unsigned z = zInner; z != zOuter; z += direction) {
                if (grid(x, y, zInner - direction)[FieldNames::VELOCITY_Z] * direction > 0.0) {
                    grid(x, y, z)[field] = grid(x, y, z - direction)[field];
                } else if (isParallelVelocity) {
                    grid(x, y, z)[field] = -grid(x, y, 2 * zInner - z - direction);
                } else {
                    grid(x, y, z)[field] = grid(x, y, 2 * zInner - z - direction);
                }
            }
        }
    }
}

template <unsigned padding>
void apply(PaddedGrid<FieldStruct, padding>& grid, const unsigned field, const unsigned face) {
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
}; // namespace Outflow
