#ifndef BASE_SOLVER
#define BASE_SOLVER

#include "../grid/base-grid.h"

class Solver {
  public:
    virtual void compute_step(Grid* grid, const double time_step) = 0;
};

#endif
