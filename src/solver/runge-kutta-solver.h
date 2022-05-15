#pragma once

#include <algorithm>
#include <optional>

#include "../boundary/boundary.h"
#include "../configuration/constants.h"
#include "../field-wrapper/fields.h"
#include "../grid/grid-functions.h"
#include "../grid/padded-grid.h"
#include "../grid/simple-grid.h"
#include "../misc/direction.h"
#include "../misc/faces.h"
#include "../misc/transformations.h"
#include "../riemann/reconstruction.h"
#include "../riemann/riemann-solver.h"
#include "../solver/base-solver.h"

typedef std::array<FieldStruct, Direction::DirMax> Changes;

template <class ProblemType>
class RungeKuttaSolver
    : public Solver<RungeKuttaSolver<ProblemType>, FieldStruct, ProblemType, GHOST_CELLS> {
  public:
    RungeKuttaSolver(PaddedGrid<FieldStruct, GHOST_CELLS>& grid,
                     const Problem<ProblemType>& problem);

    void singleStep();
    bool isFinished() const;

  private:
    SimpleGrid<FieldStruct> changeBuffer;
    SimpleGrid<FieldStruct> gridSubstepBuffer;

    double cfl;
    const unsigned rungeKuttaSteps = 2;

    void computeStep();
    void saveGrid();
    void prepareSubstep();
    void computeSubstep();
    Changes computeChanges(const unsigned x, const unsigned y, const unsigned z);
    void updateCFL(std::pair<double, double> characVelocities, const unsigned direction);
    void applyChanges(const FieldStruct& numericalValuesMinus,
                      const FieldStruct& numericalValuesPlus, const unsigned x, const unsigned y,
                      const unsigned z, const unsigned direction);
    void finaliseSubstep(const unsigned substep);
    void integrateTime(const unsigned substep);
    void advanceTime(const unsigned substep);
    void adjustTimeDelta();
    void checkErrors();
};

template <class ProblemType>
RungeKuttaSolver<ProblemType>::RungeKuttaSolver(PaddedGrid<FieldStruct, GHOST_CELLS>& grid,
                                                const Problem<ProblemType>& problem)
    : Solver<RungeKuttaSolver<ProblemType>, FieldStruct, ProblemType, GHOST_CELLS>(grid, problem),
      changeBuffer(grid.defaultValue, grid.xDim(), grid.yDim(), grid.zDim()),
      gridSubstepBuffer(grid.defaultValue, grid.xDim(), grid.yDim(), grid.zDim()) {}

template <class ProblemType> void RungeKuttaSolver<ProblemType>::singleStep() {
    this->computeStep();
    this->adjustTimeDelta();
    this->timeStep++;
}

template <class ProblemType> bool RungeKuttaSolver<ProblemType>::isFinished() const {
    return this->timeCurrent >= this->timeEnd;
}

template <class ProblemType> void RungeKuttaSolver<ProblemType>::computeStep() {
    for (unsigned substep = 0; substep < this->rungeKuttaSteps; substep++) {
        this->prepareSubstep();
        this->computeSubstep();
        this->finaliseSubstep(substep);
        this->advanceTime(substep);
    }
}

template <class ProblemType> void RungeKuttaSolver<ProblemType>::saveGrid() {
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

template <class ProblemType> void RungeKuttaSolver<ProblemType>::prepareSubstep() {
    // Number of negative temperatures is printed here, seems unnecessary

    // Start clock(s) -- time measurement omitted for now

    this->changeBuffer.clear();

    // CarbuncleFlag computation (not included by default)
}

template <class ProblemType> void RungeKuttaSolver<ProblemType>::computeSubstep() {
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

template <class ProblemType>
Changes RungeKuttaSolver<ProblemType>::computeChanges(const unsigned x, const unsigned y,
                                                      const unsigned z) {
    PerFaceValues reconstruction = Reconstruction::reconstruct(this->grid, x, y, z);

    std::array<PhysValues, Faces::FaceMax> physicalValues;
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        Transformation::reconstToConservatives(physicalValues[face], reconstruction[face],
                                               this->problem.thermal, this->problem.gamma);
        physicalValues[face].thermalPressure = Transformation::computeThermalPressure(
            reconstruction[face], this->problem.thermal, this->problem.gamma);
        RiemannSolver::computeFluxes(physicalValues[face], reconstruction[face], face);
    }

    Changes changes;
    for (unsigned dir = 0; dir < Direction::DirMax; dir++) {
        unsigned face = dir * 2;
        std::pair<double, double> characVelocities = RiemannSolver::characteristicVelocity(
            physicalValues[face], physicalValues[face + 1], reconstruction[face],
            reconstruction[face + 1], this->problem.gamma, dir);
        this->updateCFL(characVelocities, dir);
        changes[dir] = RiemannSolver::numericalFlux(characVelocities, physicalValues[face],
                                                    physicalValues[face + 1], reconstruction[face],
                                                    reconstruction[face + 1], dir);
    }
    return changes;
}

template <class ProblemType>
void RungeKuttaSolver<ProblemType>::updateCFL(const std::pair<double, double> characVelocities,
                                              const unsigned direction) {
    double maxVelocity = std::max(characVelocities.first, characVelocities.second);
    double localCFL = maxVelocity * this->problem.inverseCellSize[direction];
    this->cfl = std::max(this->cfl, localCFL);
}

template <class ProblemType>
void RungeKuttaSolver<ProblemType>::applyChanges(const FieldStruct& numericalValuesMinus,
                                                 const FieldStruct& numericalValuesPlus,
                                                 const unsigned x, const unsigned y,
                                                 const unsigned z, const unsigned direction) {
    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
        this->changeBuffer(x, y, z)[field] +=
            (numericalValuesPlus[field] - numericalValuesMinus[field]) *
            this->problem.inverseCellSize[direction];
    }
}

template <class ProblemType>
void RungeKuttaSolver<ProblemType>::finaliseSubstep(const unsigned substep) {
    this->checkErrors();

    this->problem.applySource(this->grid);
    this->checkErrors();

    Transformation::primitiveToConservative(this->grid, this->problem.thermal, this->problem.gamma);
    this->integrateTime(substep);
    Transformation::conservativeToPrimitive(this->grid, this->problem.thermal, this->problem.gamma);
    Boundary::applyAll(this->grid, this->problem);

    // Stop clock(s)
    // runtime estimation -- time measurement omitted for now

    // phystest
}

template <class ProblemType>
void RungeKuttaSolver<ProblemType>::integrateTime(const unsigned substep) {
    if (substep == 0) {
        this->saveGrid();
    }
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd(); x++) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd(); y++) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd(); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    if (substep == 0) { // if-statement should be pulled out by compiler
                        this->grid(x, y, z)[field] -=
                            this->timeDelta * this->changeBuffer(x, y, z)[field];
                    } else {
                        this->grid(x, y, z)[field] =
                            0.5 * this->gridSubstepBuffer(x, y, z)[field] +
                            0.5 * this->grid(x, y, z)[field] -
                            0.5 * this->timeDelta * this->changeBuffer(x, y, z)[field];
                    }
                }
            }
        }
    }
}

template <class ProblemType>
void RungeKuttaSolver<ProblemType>::advanceTime(const unsigned substep) {
    if (substep == 0) {
        this->timeCurrent += this->timeDelta;
    }
}

template <class ProblemType> void RungeKuttaSolver<ProblemType>::adjustTimeDelta() {
    this->timeDelta = this->problem.cflThreshold / this->cfl;
}

template <class ProblemType> void RungeKuttaSolver<ProblemType>::checkErrors() {
    if (GridFunctions::checkNaN(this->grid)) {
        std::cerr << "Encountered NaN" << std::endl;
    }
}
