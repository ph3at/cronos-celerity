#pragma once

#include <omp.h>

#include "../boundary/boundary.h"
#include "../data-types/direction.h"
#include "../data-types/faces.h"
#include "../grid/grid-functions.h"
#include "../grid/padded-grid.h"
#include "../grid/simple-grid.h"
#include "../riemann/reconstruction.h"
#include "../riemann/riemann-solver.h"
#include "../solver/base-solver.h"
#include "../transformation/transformations.h"

template <class Fields> using Changes = std::array<Fields, Direction::DirMax>;

template <class ProblemType, class Fields, unsigned padding>
class RungeKuttaSolver
    : public Solver<RungeKuttaSolver<ProblemType, Fields, padding>, Fields, ProblemType, padding> {
  public:
    RungeKuttaSolver(const ProblemType& problem, const unsigned rungeKuttaSteps = 2,
                     const bool doOutput = false);
    RungeKuttaSolver(const PaddedGrid<Fields, padding>& grid, const ProblemType& problem,
                     const unsigned rungeKuttaSteps = 2, const bool doOutput = false);

    void init();
    void singleStep();
    void adjustConfig();
    void finaliseResult(){};

  private:
    SimpleGrid<Changes<Fields>> changeBuffer;
    SimpleGrid<Fields> gridSubstepBuffer;

    double cfl;
    const unsigned rungeKuttaSteps;

    void saveGrid();
    void prepareSubstep();
    void computeSubstep();
    Changes<Fields> computeChanges(const unsigned x, const unsigned y, const unsigned z);
    void updateCFL(std::pair<double, double> characVelocities, const unsigned direction);
    void finaliseSubstep(const unsigned substep);
    void integrateTime(const unsigned substep);
    void checkErrors();
};

template <class ProblemType, class Fields, unsigned padding>
RungeKuttaSolver<ProblemType, Fields, padding>::RungeKuttaSolver(const ProblemType& problem,
                                                                 const unsigned rungeKuttaSteps,
                                                                 const bool doOutput)
    : Solver<RungeKuttaSolver<ProblemType, Fields, padding>, Fields, ProblemType, padding>(
          problem, doOutput),
      changeBuffer({}, problem.numberCells[Direction::DirX] + 2 * padding,
                   problem.numberCells[Direction::DirY] + 2 * padding,
                   problem.numberCells[Direction::DirZ] + 2 * padding),
      gridSubstepBuffer({}, problem.numberCells[Direction::DirX] + 2 * padding,
                        problem.numberCells[Direction::DirY] + 2 * padding,
                        problem.numberCells[Direction::DirZ] + 2 * padding),
      rungeKuttaSteps(rungeKuttaSteps) {}

template <class ProblemType, class Fields, unsigned padding>
RungeKuttaSolver<ProblemType, Fields, padding>::RungeKuttaSolver(
    const PaddedGrid<Fields, padding>& grid, const ProblemType& problem,
    const unsigned rungeKuttaSteps, const bool doOutput)
    : Solver<RungeKuttaSolver<ProblemType, Fields, padding>, Fields, ProblemType, padding>(
          grid, problem, doOutput),
      changeBuffer({}, grid.xDim(), grid.yDim(), grid.zDim()),
      gridSubstepBuffer(grid.defaultValue, grid.xDim(), grid.yDim(), grid.zDim()),
      rungeKuttaSteps(rungeKuttaSteps) {}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::init() {
    this->problem.initialiseGrid(this->grid);
    Boundary::applyAll(this->grid, this->problem);
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::singleStep() {
    this->cfl = 0.0;
    for (unsigned substep = 0; substep < this->rungeKuttaSteps; substep++) {
        this->prepareSubstep();
        this->computeSubstep();
        this->finaliseSubstep(substep);
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::saveGrid() {
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd(); x++) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd(); y++) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd(); z++) {
                for (unsigned field = 0; field < Fields().size(); field++) {
                    this->gridSubstepBuffer(x, y, z)[field] = this->grid(x, y, z)[field];
                }
            }
        }
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::prepareSubstep() {
    // Number of negative temperatures is printed here, seems unnecessary

    // Start clock(s) -- time measurement omitted for now

    this->changeBuffer.clear();

    // CarbuncleFlag computation (not included by default)
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::computeSubstep() {
#pragma omp parallel for
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd() + 1; x++) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd() + 1; y++) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd() + 1; z++) {
                this->changeBuffer(x, y, z) = this->computeChanges(x, y, z);
            }
        }
    }
}

template <class ProblemType, class Fields, unsigned padding>
Changes<Fields> RungeKuttaSolver<ProblemType, Fields, padding>::computeChanges(const unsigned x,
                                                                               const unsigned y,
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

    Changes<Fields> changes;
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

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::updateCFL(
    const std::pair<double, double> characVelocities, const unsigned direction) {
    double maxVelocity = std::max(characVelocities.first, characVelocities.second);
    double localCFL = maxVelocity * this->problem.inverseCellSize[direction];
    this->cfl = std::max(this->cfl, localCFL);
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::finaliseSubstep(const unsigned substep) {
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

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::integrateTime(const unsigned substep) {
    if (substep == 0 && this->rungeKuttaSteps > 1) {
        this->saveGrid();
    }
#pragma omp parallel for
    for (unsigned x = this->grid.xStart(); x < this->grid.xEnd(); x++) {
        for (unsigned y = this->grid.yStart(); y < this->grid.yEnd(); y++) {
            for (unsigned z = this->grid.zStart(); z < this->grid.zEnd(); z++) {
                for (unsigned field = 0; field < Fields().size(); field++) {
                    const double changeX =
                        (this->changeBuffer(x + 1, y, z)[Direction::DirX][field] -
                         this->changeBuffer(x, y, z)[Direction::DirX][field]) *
                        this->problem.inverseCellSize[Direction::DirX];
                    const double changeY =
                        (this->changeBuffer(x, y + 1, z)[Direction::DirY][field] -
                         this->changeBuffer(x, y, z)[Direction::DirY][field]) *
                        this->problem.inverseCellSize[Direction::DirY];
                    const double changeZ =
                        (this->changeBuffer(x, y, z + 1)[Direction::DirZ][field] -
                         this->changeBuffer(x, y, z)[Direction::DirZ][field]) *
                        this->problem.inverseCellSize[Direction::DirZ];
                    const double change = changeX + changeY + changeZ;
                    if (substep == 0) { // if-statement should be pulled out by compiler
                        this->grid(x, y, z)[field] -= this->timeDelta * change;
                    } else {
                        this->grid(x, y, z)[field] = 0.5 * this->gridSubstepBuffer(x, y, z)[field] +
                                                     0.5 * this->grid(x, y, z)[field] -
                                                     0.5 * this->timeDelta * change;
                    }
                }
            }
        }
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::adjustConfig() {
    this->timeDelta = this->problem.cflThreshold / this->cfl;
    if (this->problem.preciseEnd) {
        this->timeDelta = std::min(this->timeDelta, this->timeEnd - this->timeCurrent);
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSolver<ProblemType, Fields, padding>::checkErrors() {
    if (GridFunctions::checkNaN(this->grid)) {
        std::cerr << "Encountered NaN" << std::endl;
    }
}
