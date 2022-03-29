#pragma once

#include "../grid/base-grid.h"

template <class T> class Solver {
  public:
    virtual void solve() = 0;

  protected:
    Solver(Grid<T>& grid, const double timeDelta, const double timeStart, const double timeEnd);

    Grid<T>& grid;
    double timeDelta;
    const double timeStart = 0.0;
    const double timeEnd;
};

template <class T>
Solver<T>::Solver(Grid<T>& grid, const double timeDelta, const double timeStart,
                  const double timeEnd)
    : grid(grid), timeStart(timeStart), timeEnd(timeEnd) {
    this->timeDelta = timeDelta;
}
