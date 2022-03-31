#include "runge-kutta-solver.h"

RungeKuttaSolver::RungeKuttaSolver(PaddedGrid<FieldStruct, GHOST_CELLS>& grid,
                                   const double timeDelta, const double timeStart,
                                   const double timeEnd, const unsigned rungeKuttaSteps)
    : Solver<FieldStruct, GHOST_CELLS>(grid, timeDelta, timeStart, timeEnd),
      changeBuffer(grid.defaultValue, grid.xDim(), grid.yDim(), grid.zDim()),
      rungeKuttaSteps(rungeKuttaSteps) {
    this->timeCurrent = timeStart;
    this->timeStep = 0;
};

void RungeKuttaSolver::solve() {
    while (!this->isFinished()) {
        this->computeStep();
        this->adjustTimeDelta();
        this->timeStep++;
    }
}

bool RungeKuttaSolver::isFinished() const {
    return this->timeCurrent >= this->timeEnd;
}

void RungeKuttaSolver::computeStep() {
    for (unsigned substep = 0; substep < this->rungeKuttaSteps; substep++) {
        this->prepareSubstep();
        this->computeSubstep();
        this->finaliseSubstep();
    }
}

void RungeKuttaSolver::prepareSubstep() {
    this->changeBuffer.clear();
}

void RungeKuttaSolver::computeSubstep() {
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd(); x++) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd(); y++) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd(); z++) {
                // TODO: Find more descriptive name than numVals
                const PhysValues numVals = this->computeChanges(x, y, z);

                const PhysValues numValsX = this->computeChanges(x + 1, y, z);
                this->applyChanges(Direction::DirX, x, y, z, numVals, numValsX);

                const PhysValues numValsY = this->computeChanges(x, y + 1, z);
                this->applyChanges(Direction::DirY, x, y, z, numVals, numValsY);

                const PhysValues numValsZ = this->computeChanges(x, y, z + 1);
                this->applyChanges(Direction::DirZ, x, y, z, numVals, numValsZ);
            }
        }
    }
}

PhysValues RungeKuttaSolver::computeChanges(const unsigned x, const unsigned y, const unsigned z) {
    PhysValues reconstruction = this->reconstruct(x, y, z);
    return {};
}

PhysValues RungeKuttaSolver::reconstruct(const unsigned x, const unsigned y, const unsigned z) {
    return {};
}

void RungeKuttaSolver::applyChanges(const Direction direction, const unsigned x, const unsigned y,
                                    const unsigned z, const PhysValues numVals,
                                    const PhysValues numValsX) {}

void RungeKuttaSolver::finaliseSubstep() {
    this->timeCurrent += this->timeDelta;
}

void RungeKuttaSolver::adjustTimeDelta() {}
