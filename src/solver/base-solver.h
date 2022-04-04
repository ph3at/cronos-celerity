#pragma once

#include "../grid/padded-grid.h"
#include "../parameters/problem.h"

template <class T, unsigned ghostCells> class Solver {
  public:
    virtual void solve() = 0;

  protected:
    Solver(PaddedGrid<T, ghostCells>& grid, const Problem& problem);

    PaddedGrid<T, ghostCells>& grid;
    const Problem& problem;
    double timeDelta;
};

template <class T, unsigned ghostCells>
Solver<T, ghostCells>::Solver(PaddedGrid<T, ghostCells>& grid, const Problem& problem)
    : grid(grid), problem(problem) {
    this->timeDelta = problem.timeDelta;
}
