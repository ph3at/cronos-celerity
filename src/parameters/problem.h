#pragma once

#include <array>

#include "../misc/direction.h"

typedef struct problem {
    const bool thermal;
    const double timeDelta;
    const double timeStart;
    const double timeEnd;
    const double gamma;
    const std::array<std::size_t, Direction::DirMax> numberCells;
    const std::array<double, Direction::DirMax> cellSize;
    const std::array<double, Direction::DirMax> inverseCellSize;

    static struct problem initialiseTestProblem();
} Problem;
