#pragma once

#include <vector>

#include "../boundary/boundary.h"
#include "../data-types/direction.h"
#include "../data-types/faces.h"
#include "../grid/grid-functions.h"
#include "../grid/padded-grid.h"
#include "../riemann/reconstruction.h"
#include "../riemann/riemann-solver.h"
#include "../transformation/transformations.h"

template <class Fields> using Changes = std::array<Fields, Direction::DirMax>;

template <class ProblemType, class Fields, unsigned padding> class RungeKuttaSyclSolver {
  public:
    RungeKuttaSyclSolver(const ProblemType& problem, const unsigned rungeKuttaSteps = 2);
    RungeKuttaSyclSolver(const PaddedGrid<Fields, padding>& grid, const ProblemType& problem,
                         const unsigned rungeKuttaSteps = 2);

    void initialise();
    void step();
    void adjust();
    void finaliseResult(){};

    PaddedGrid<Fields, padding> grid;

  private:
    const std::size_t m_sizeX;
    const std::size_t m_sizeY;
    const std::size_t m_sizeZ;

    std::vector<Changes<Fields>> changeBuffer;
    std::vector<Fields> gridSubstepBuffer;

    double cfl;
    const unsigned rungeKuttaSteps;

    const ProblemType& problem;

    double timeDelta;
    double timeCurrent;
    const double timeEnd;

    void saveGrid();
    void prepareSubstep();
    void computeSubstep();
    Changes<Fields> computeChanges(const unsigned x, const unsigned y, const unsigned z);
    void updateCFL(std::pair<double, double> characVelocities, const unsigned direction);
    void finaliseSubstep(const unsigned substep);
    void integrateTime(const unsigned substep);
    void checkErrors();

    std::size_t idx3d(const std::size_t x, const std::size_t y,
                      const std::size_t z) const noexcept {
        return x * m_sizeY * m_sizeZ + y * m_sizeZ + z;
    }
};

template <class ProblemType, class Fields, unsigned padding>
RungeKuttaSyclSolver<ProblemType, Fields, padding>::RungeKuttaSyclSolver(
    const ProblemType& problem, const unsigned rungeKuttaSteps)
    : grid({}, problem.numberCells[Direction::DirX], problem.numberCells[Direction::DirY],
           problem.numberCells[Direction::DirZ]),
      m_sizeX(problem.numberCells[Direction::DirX] + 2 * padding),
      m_sizeY(problem.numberCells[Direction::DirY] + 2 * padding),
      m_sizeZ(problem.numberCells[Direction::DirZ] + 2 * padding),
      changeBuffer(m_sizeX * m_sizeY * m_sizeZ), gridSubstepBuffer(m_sizeX * m_sizeY * m_sizeZ),
      rungeKuttaSteps(rungeKuttaSteps), problem(problem), timeDelta(problem.timeDelta),
      timeCurrent(problem.timeStart), timeEnd(problem.timeEnd) {}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::initialise() {
    problem.initialiseGrid(grid);
    Boundary::applyAll(grid, problem);
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::step() {
    cfl = 0.0;
    for (unsigned substep = 0; substep < rungeKuttaSteps; substep++) {
        prepareSubstep();
        computeSubstep();
        finaliseSubstep(substep);
    }
    timeCurrent += timeDelta;
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::saveGrid() {
    for (unsigned x = grid.xStart(); x < grid.xEnd(); x++) {
        for (unsigned y = grid.yStart(); y < grid.yEnd(); y++) {
            for (unsigned z = grid.zStart(); z < grid.zEnd(); z++) {
                for (unsigned field = 0; field < Fields().size(); field++) {
                    gridSubstepBuffer[idx3d(x, y, z)][field] = grid(x, y, z)[field];
                }
            }
        }
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::prepareSubstep() {
    // Number of negative temperatures is printed here, seems unnecessary

    // Start clock(s) -- time measurement omitted for now

    changeBuffer.clear();

    // CarbuncleFlag computation (not included by default)
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::computeSubstep() {
    for (unsigned x = grid.xStart(); x < grid.xEnd() + 1; x++) {
        for (unsigned y = grid.yStart(); y < grid.yEnd() + 1; y++) {
            for (unsigned z = grid.zStart(); z < grid.zEnd() + 1; z++) {
                changeBuffer[idx3d(x, y, z)] = computeChanges(x, y, z);
            }
        }
    }
}

template <class ProblemType, class Fields, unsigned padding>
Changes<Fields> RungeKuttaSyclSolver<ProblemType, Fields, padding>::computeChanges(
    const unsigned x, const unsigned y, const unsigned z) {
    PerFaceValues reconstruction = Reconstruction::reconstruct(grid, x, y, z);

    std::array<PhysValues, Faces::FaceMax> physicalValues;
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        Transformation::reconstToConservatives(physicalValues[face], reconstruction[face],
                                               problem.thermal, problem.gamma);
        physicalValues[face].thermalPressure = Transformation::computeThermalPressure(
            reconstruction[face], problem.thermal, problem.gamma);
        RiemannSolver::computeFluxes(physicalValues[face], reconstruction[face], face);
    }

    Changes<Fields> changes;
    for (unsigned dir = 0; dir < Direction::DirMax; dir++) {
        unsigned face = dir * 2;
        std::pair<double, double> characVelocities = RiemannSolver::characteristicVelocity(
            physicalValues[face], physicalValues[face + 1], reconstruction[face],
            reconstruction[face + 1], problem.gamma, dir);
        updateCFL(characVelocities, dir);
        changes[dir] = RiemannSolver::numericalFlux(characVelocities, physicalValues[face],
                                                    physicalValues[face + 1], reconstruction[face],
                                                    reconstruction[face + 1], dir);
    }
    return changes;
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::updateCFL(
    const std::pair<double, double> characVelocities, const unsigned direction) {
    double maxVelocity = std::max(characVelocities.first, characVelocities.second);
    double localCFL = maxVelocity * problem.inverseCellSize[direction];
    cfl = std::max(cfl, localCFL);
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::finaliseSubstep(const unsigned substep) {
    checkErrors();

    problem.applySource(grid);
    checkErrors();

    Transformation::primitiveToConservative(grid, problem.thermal, problem.gamma);
    integrateTime(substep);
    Transformation::conservativeToPrimitive(grid, problem.thermal, problem.gamma);

    Boundary::applyAll(grid, problem);

    // Stop clock(s)
    // runtime estimation -- time measurement omitted for now

    // phystest
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::integrateTime(const unsigned substep) {
    if (substep == 0 && rungeKuttaSteps > 1) {
        saveGrid();
    }
    for (unsigned x = grid.xStart(); x < grid.xEnd(); x++) {
        for (unsigned y = grid.yStart(); y < grid.yEnd(); y++) {
            for (unsigned z = grid.zStart(); z < grid.zEnd(); z++) {
                for (unsigned field = 0; field < Fields().size(); field++) {
                    const double changeX =
                        (changeBuffer[idx3d(x + 1, y, z)][Direction::DirX][field] -
                         changeBuffer[idx3d(x, y, z)][Direction::DirX][field]) *
                        problem.inverseCellSize[Direction::DirX];
                    const double changeY =
                        (changeBuffer[idx3d(x, y + 1, z)][Direction::DirY][field] -
                         changeBuffer[idx3d(x, y, z)][Direction::DirY][field]) *
                        problem.inverseCellSize[Direction::DirY];
                    const double changeZ =
                        (changeBuffer[idx3d(x, y, z + 1)][Direction::DirZ][field] -
                         changeBuffer[idx3d(x, y, z)][Direction::DirZ][field]) *
                        problem.inverseCellSize[Direction::DirZ];
                    const double change = changeX + changeY + changeZ;
                    if (substep == 0) { // if-statement should be pulled out by compiler
                        grid(x, y, z)[field] -= timeDelta * change;
                    } else {
                        grid(x, y, z)[field] = 0.5 * gridSubstepBuffer[idx3d(x, y, z)][field] +
                                               0.5 * grid(x, y, z)[field] -
                                               0.5 * timeDelta * change;
                    }
                }
            }
        }
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::adjust() {
    timeDelta = problem.cflThreshold / cfl;
    if (problem.preciseEnd) {
        timeDelta = std::min(timeDelta, timeEnd - timeCurrent);
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::checkErrors() {
    if (GridFunctions::checkNaN(grid)) {
        std::cerr << "Encountered NaN" << std::endl;
    }
}
