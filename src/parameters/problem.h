#pragma once

#include <array>

#include "../boundary/boundary-types.h"
#include "../misc/direction.h"
#include "../misc/faces.h"

typedef struct problem {
    const bool thermal;
    const double timeDelta;
    const double timeStart;
    const double timeEnd;
    const double gamma;
    const std::array<std::size_t, Direction::DirMax> numberCells;
    const std::array<double, Direction::DirMax> cellSize;
    const std::array<double, Direction::DirMax> inverseCellSize;

    const std::array<BoundaryType, Faces::FaceMax> boundaryTypes;

    static struct problem initialiseTestProblem();
} Problem;
