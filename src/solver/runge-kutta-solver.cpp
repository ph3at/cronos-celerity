#include "runge-kutta-solver.h"
#include "../misc/transformations.h"
#include "../reconstruction/reconstruction.h"
#include "../riemann/riemann-solver.h"

RungeKuttaSolver::RungeKuttaSolver(PaddedGrid<FieldStruct, GHOST_CELLS>& grid,
                                   const Problem& problem, const unsigned rungeKuttaSteps)
    : Solver<FieldStruct, GHOST_CELLS>(grid, problem),
      changeBuffer(grid.defaultValue, grid.xDim(), grid.yDim(), grid.zDim()),
      rungeKuttaSteps(rungeKuttaSteps) {
    this->timeCurrent = problem.timeStart;
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
    return this->timeCurrent >= this->problem.timeEnd;
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
                const Changes numVals = this->computeChanges(x, y, z);

                const Changes numValsX = this->computeChanges(x + 1, y, z);
                this->applyChanges(Direction::DirX, x, y, z, numVals, numValsX);

                const Changes numValsY = this->computeChanges(x, y + 1, z);
                this->applyChanges(Direction::DirY, x, y, z, numVals, numValsY);

                const Changes numValsZ = this->computeChanges(x, y, z + 1);
                this->applyChanges(Direction::DirZ, x, y, z, numVals, numValsZ);
            }
        }
    }
}

Changes RungeKuttaSolver::computeChanges(const unsigned x, const unsigned y, const unsigned z) {
    PerFaceValues reconstruction = Reconstruction::reconstruct(this->grid, x, y, z);

    std::array<PhysValues, Faces::FaceMax> physicalValues;
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        Transformation::reconstToConservatives(physicalValues[face], reconstruction[face],
                                               this->problem);
        physicalValues[face].thermalPressure =
            Transformation::computeThermalPressure(reconstruction[face], this->problem);
        RiemannSolver::computeFluxes(physicalValues[face], reconstruction[face], face);
    }
    Changes changes;
    for (unsigned dir = 0; dir < Direction::DirMax; dir++) {
        unsigned face = dir * 2;
        std::pair<double, double> characVelocities = RiemannSolver::characteristicVelocity(
            physicalValues[face], physicalValues[face + 1], reconstruction[face],
            reconstruction[face + 1], this->problem, dir);
        this->updateCFL(characVelocities);
        changes[dir] = RiemannSolver::numericalFlux(characVelocities, physicalValues[face],
                                                    physicalValues[face + 1], dir);
    }
    return changes;
}

void RungeKuttaSolver::updateCFL(const std::pair<double, double> characVelocities) {
    double maxVelocity = std::max(characVelocities.first, characVelocities.second);
    double localCFL = maxVelocity * 0.1; // TODO: 0.1 -> inverseDX[dir]
    this->cfl = std::max(this->cfl, localCFL);
}

void RungeKuttaSolver::applyChanges(const Direction direction, const unsigned x, const unsigned y,
                                    const unsigned z, const Changes numVals,
                                    const Changes numValsX) {}

void RungeKuttaSolver::finaliseSubstep() {
    this->timeCurrent += this->timeDelta;
}

void RungeKuttaSolver::adjustTimeDelta() {}
