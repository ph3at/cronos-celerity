#pragma once

#include "base-solver.h"

class SimpleSolver : public Solver<double> {
  public:
    SimpleSolver(Grid<double>& grid, const double timeDelta, const double timeEnd)
        : Solver(grid, timeDelta, timeEnd){};

    void solve();
    void computeStep();

  private:
    SimpleSolver();
};
