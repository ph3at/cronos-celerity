#pragma once

#include <algorithm>

#include "../field-wrapper/fields.h"
#include "../grid/padded-grid.h"
#include "../grid/simple-grid.h"
#include "../misc/direction.h"
#include "../misc/faces.h"
#include "../parameters/constants.h"
#include "../solver/base-solver.h"

typedef std::array<FieldStruct, Direction::DirMax> Changes;

class RungeKuttaSolver : public Solver<FieldStruct, GHOST_CELLS> {
  public:
    RungeKuttaSolver(PaddedGrid<FieldStruct, GHOST_CELLS>& grid, const Problem& problem,
                     unsigned rungeKuttaSteps);

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
    Changes computeChanges(const unsigned x, const unsigned y, const unsigned z);
    void updateCFL(std::pair<double, double> characVelocities, const unsigned direction);
    void applyChanges(const FieldStruct& numericalValuesMinus,
                      const FieldStruct& numericalValuesPlus, const unsigned x, const unsigned y,
                      const unsigned z, const unsigned direction);
    void finaliseSubstep();
    void adjustTimeDelta();
};
