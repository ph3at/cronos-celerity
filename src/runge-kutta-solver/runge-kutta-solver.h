#pragma once

#include <algorithm>

#include "../grid/padded-grid.h"
#include "../grid/simple-grid.h"
#include "../misc/direction.h"
#include "../misc/phys-fields.h"
#include "../solver/base-solver.h"

typedef std::vector<double> PhysValues;

class RungeKuttaSolver : public Solver<FieldStruct> {
  public:
    RungeKuttaSolver(PaddedGrid<FieldStruct>& grid, const double timeDelta, const double timeStart,
                     const double timeEnd, unsigned rungeKuttaSteps);

    void solve();

  private:
    SimpleGrid<FieldStruct> changeBuffer;

    double timeCurrent;
    unsigned timeStep;
    double cfl;
    const unsigned rungeKuttaSteps = 2;

    void computeStep();
    bool isFinished() const;
    void prepareSubstep();
    void computeSubstep();
    PhysValues computeChanges(const unsigned x, const unsigned y, const unsigned z);
    PhysValues reconstruct(const unsigned x, const unsigned y, const unsigned z);
    void applyChanges(const Direction direction, const unsigned x, const unsigned y,
                      const unsigned z, const PhysValues numVals, const PhysValues numValsX);
    void finaliseSubstep();
    void adjustTimeDelta();
};
