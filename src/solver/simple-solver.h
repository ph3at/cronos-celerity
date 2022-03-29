#pragma once

#include "base-solver.h"

class SimpleSolver : public Solver<double> {
  public:
    SimpleSolver(Grid<double>& grid, const double timeDelta, const double timeEnd)
        : Solver(grid, timeDelta, 0.0, timeEnd){};

    void solve();
    void computeStep();

  private:
    SimpleSolver();
};
