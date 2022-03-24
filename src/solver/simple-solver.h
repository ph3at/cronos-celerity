#pragma once

#include "base-solver.h"

class SimpleSolver : public Solver<double> {
  public:
    void computeStep(Grid<double>& grid, const double timeStep);
};
