#pragma once

#include "../field-wrapper/phys-fields.h"
#include "../grid/padded-grid.h"
#include "../misc/faces.h"

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
                grid(x, y, z)[field] = grid(x - direction, y, z)[field];
            }
        }
    }
}

template <unsigned padding>
void Extrapolate::applyY(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                         const bool isLeft) {
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
                grid(x, y, z)[field] = grid(x, y - direction, z)[field];
            }
        }
    }
}

template <unsigned padding>
void Extrapolate::applyZ(PaddedGrid<FieldStruct, padding>& grid, const unsigned field,
                         const bool isLeft) {
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

