#ifndef SIMPLE_SOLVER
#define SIMPLE_SOLVER

#include "base-solver.h"

class SimpleSolver : public Solver {
  public:
    void computeStep(Grid& grid, const double timeStep);
};

#endif
