#include "runge-kutta-solver.h"
#include "../boundary/boundary.h"
#include "../grid/grid-functions.h"
#include "../misc/transformations.h"
#include "../reconstruction/reconstruction.h"
#include "../riemann/riemann-solver.h"

RungeKuttaSolver::RungeKuttaSolver(PaddedGrid<FieldStruct, GHOST_CELLS>& grid,
                                   const Problem& problem)
    : Solver<FieldStruct, GHOST_CELLS>(grid, problem),
      changeBuffer(grid.defaultValue, grid.xDim(), grid.yDim(), grid.zDim()),
      gridSubstepBuffer(grid.defaultValue, grid.xDim(), grid.yDim(), grid.zDim()) {
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
    this->saveGrid();
    for (unsigned substep = 0; substep < this->rungeKuttaSteps; substep++) {
        this->prepareSubstep();
        this->computeSubstep();
        this->finaliseSubstep(substep);
        this->advanceTime(substep);
    }
}

void RungeKuttaSolver::saveGrid() {
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd(); x++) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd(); y++) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd(); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    this->gridSubstepBuffer(x, y, z)[field] = this->grid(x, y, z)[field];
                }
            }
        }
    }
}

void RungeKuttaSolver::prepareSubstep() {
    // Number of negative temperatures is printed here, seems unnecessary

    // Start clock(s) -- time measurement omitted for now

    this->changeBuffer.clear();

    /* This part seems to never be used, so I will not implement it for now.
     *
     * Transformation::primitiveToConservative(this->grid, this->problem);
     * Save old variables
     * Transformation::conservativeToPrimitive
     * grid.checkNaN() */

    // CarbuncleFlag computation (not included by default)
}

void RungeKuttaSolver::computeSubstep() {
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd(); x++) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd(); y++) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd(); z++) {
                const Changes numericalFluxesCenter = this->computeChanges(x, y, z);

                const Changes numericalFluxesX = this->computeChanges(x + 1, y, z);
                this->applyChanges(numericalFluxesCenter[Direction::DirX],
                                   numericalFluxesX[Direction::DirX], x, y, z, Direction::DirX);

                const Changes numericalFluxesY = this->computeChanges(x, y + 1, z);
                this->applyChanges(numericalFluxesCenter[Direction::DirY],
                                   numericalFluxesY[Direction::DirY], x, y, z, Direction::DirY);

                const Changes numericalFluxesZ = this->computeChanges(x, y, z + 1);
                this->applyChanges(numericalFluxesCenter[Direction::DirZ],
                                   numericalFluxesZ[Direction::DirZ], x, y, z, Direction::DirZ);
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
        this->updateCFL(characVelocities, dir);
        changes[dir] = RiemannSolver::numericalFlux(characVelocities, physicalValues[face],
                                                    physicalValues[face + 1], reconstruction[face],
                                                    reconstruction[face + 1], dir);
    }
    return changes;
}

void RungeKuttaSolver::updateCFL(const std::pair<double, double> characVelocities,
                                 const unsigned direction) {
    double maxVelocity = std::max(characVelocities.first, characVelocities.second);
    double localCFL = maxVelocity * this->problem.inverseCellSize[direction];
    this->cfl = std::max(this->cfl, localCFL);
}

void RungeKuttaSolver::applyChanges(const FieldStruct& numericalValuesMinus,
                                    const FieldStruct& numericalValuesPlus, const unsigned x,
                                    const unsigned y, const unsigned z, const unsigned direction) {
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        this->changeBuffer(x, y, z)[field] +=
            (numericalValuesPlus[field] - numericalValuesMinus[field]) *
            this->problem.inverseCellSize[direction];
    }
}

void RungeKuttaSolver::finaliseSubstep(const unsigned substep) {
    if (GridFunctions::checkNaN(this->grid)) {
        std::cerr << "Encountered NaN" << std::endl;
    }
    Transformation::primitiveToConservative(this->grid, this->problem);
    this->integrateTime(substep);
    Transformation::conservativeToPrimitive(this->grid, this->problem);
    Boundary::applyAll(this->grid, this->problem);

    // Stop clock(s)
    // runtime estimation -- time measurement omitted for now

    // phystest
}

void RungeKuttaSolver::integrateTime(const unsigned substep) {
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd(); x++) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd(); y++) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd(); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    if (substep == 0) { // if-statement should be pulled out by compiler
                        this->grid(x, y, z)[field] -=
                            this->timeDelta * changeBuffer(x, y, z)[field];
                    } else {
                        this->grid(x, y, z)[field] =
                            0.5 * this->gridSubstepBuffer(x, y, z)[field] +
                            0.5 * this->grid(x, y, z)[field] -
                            0.5 * this->timeDelta * changeBuffer(x, y, z)[field];
                    }
                }
            }
        }
    }
}

void RungeKuttaSolver::advanceTime(const unsigned substep) {
    if (substep == 0) {
        this->timeCurrent += this->timeDelta;
    }
}

void RungeKuttaSolver::adjustTimeDelta() {}
