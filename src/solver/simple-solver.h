#ifndef SIMPLE_SOLVER
#define SIMPLE_SOLVER

#include "base-solver.h"

class SimpleSolver : public Solver {
  public:
    void compute_step(Grid& grid, const double time_step);
};

#endif
