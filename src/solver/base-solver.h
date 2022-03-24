#ifndef BASE_SOLVER
#define BASE_SOLVER

#include "../grid/base-grid.h"

class Solver {
  public:
    virtual void computeStep(Grid& grid, const double timeStep) = 0;
};

#endif
