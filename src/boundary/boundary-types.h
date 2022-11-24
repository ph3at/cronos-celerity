#pragma once

enum BoundaryType { EMPTY, EXTRAPOLATE, OUTFLOW, USER };

namespace Boundary {

enum class Axis : int {
    X = 0,
    Y = 1,
    Z = 2,
};

enum class Dir {
    LEFT,
    RIGHT,
};

} // namespace Boundary
