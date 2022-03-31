#pragma once

#include "../grid/padded-grid.h"

template <class T, unsigned ghostCells> class Solver {
  public:
    virtual void solve() = 0;

  protected:
    Solver(PaddedGrid<T, ghostCells>& grid, const double timeDelta, const double timeStart,
           const double timeEnd);

    PaddedGrid<T, ghostCells>& grid;
    double timeDelta;
    const double timeStart = 0.0;
    const double timeEnd;
};

template <class T, unsigned ghostCells>
Solver<T, ghostCells>::Solver(PaddedGrid<T, ghostCells>& grid, const double timeDelta,
                              const double timeStart, const double timeEnd)
    : grid(grid), timeStart(timeStart), timeEnd(timeEnd) {
    this->timeDelta = timeDelta;
}
