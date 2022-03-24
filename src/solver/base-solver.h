#pragma once

#include "../grid/base-grid.h"

template <class T> class Solver {
  public:
    virtual void computeStep(Grid<T>& grid, const double timeStep) = 0;
};
