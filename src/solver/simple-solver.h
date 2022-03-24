#pragma once

#include "base-solver.h"

class SimpleSolver : public Solver {
  public:
    void computeStep(Grid& grid, const double timeStep);
};
