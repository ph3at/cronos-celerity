#pragma once

#include "../configuration/problem.h"
#include "../grid/padded-grid.h"

template <class Fields, class ProblemType, unsigned ghostCells> class Solver {
  public:
    virtual void solve() = 0;

  protected:
    Solver(PaddedGrid<Fields, ghostCells>& grid, const Problem<ProblemType>& problem);

    PaddedGrid<Fields, ghostCells>& grid;
    const Problem<ProblemType>& problem;
    double timeDelta;
};

template <class Fields, class ProblemType, unsigned ghostCells>
Solver<Fields, ProblemType, ghostCells>::Solver(PaddedGrid<Fields, ghostCells>& grid,
                                                const Problem<ProblemType>& problem)
    : grid(grid), problem(problem) {
    this->timeDelta = problem.timeDelta;
}
