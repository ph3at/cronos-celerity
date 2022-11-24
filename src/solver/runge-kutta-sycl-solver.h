#pragma once

#include <algorithm>
#include <tuple>
#include <vector>

#include <CL/sycl.hpp>

#include "../boundary/boundary.h"
#include "../data-types/direction.h"
#include "../data-types/faces.h"
#include "../grid/grid-functions.h"
#include "../grid/padded-grid.h"
#include "../grid/utils.h"
#include "../riemann/reconstruction.h"
#include "../riemann/riemann-solver.h"
#include "../transformation/transformations.h"

template <class Fields> using Changes = std::array<Fields, Direction::DirMax>;

template <class ProblemType, class Fields, unsigned padding> class RungeKuttaSyclSolver {
  public:
    RungeKuttaSyclSolver(const ProblemType& problem, const unsigned rungeKuttaSteps = 2);

    void initialise();
    void step();
    void adjust();
    void finaliseResult(){};

    PaddedGrid<Fields, padding> grid() const {
        auto paddedGrid =
            PaddedGrid<Fields, padding>({}, m_sizeX - 2 * padding, m_sizeY - 2 * padding, m_sizeZ - 2 * padding);
        for (std::size_t x = 0; x < m_sizeX; ++x) {
            for (std::size_t y = 0; y < m_sizeY; ++y) {
                for (std::size_t z = 0; z < m_sizeZ; ++z) {
                    paddedGrid(x, y, z) = m_grid[idx3d(x, y, z)];
                }
            }
        }
        return paddedGrid;
    }

  private:
    cl::sycl::queue m_queue;

    const std::size_t m_sizeX;
    const std::size_t m_sizeY;
    const std::size_t m_sizeZ;
    const grid::utils::dimensions m_dims;

    std::vector<Fields> m_grid;
    std::vector<Changes<Fields>> changeBuffer;
    std::vector<Fields> gridSubstepBuffer;
    std::vector<double> m_cflBuffer;

    double m_cfl;
    const unsigned rungeKuttaSteps;

    const ProblemType& problem;

    double timeDelta;
    double timeCurrent;
    const double timeEnd;

    void prepareSubstep();
    void computeSubstep();
    std::pair<Changes<Fields>, double> computeChanges(const cl::sycl::id<3>& id) const;
    double calcCFL(std::pair<double, double> characVelocities, const unsigned direction) const;
    double reduceCFL(const double initialCFL, const std::vector<double>& cflBuffer) const;
    void finaliseSubstep(const unsigned substep);
    void integrateTime(const unsigned substep);
    void checkErrors();

    std::size_t idx3d(const std::size_t x, const std::size_t y, const std::size_t z) const noexcept {
        return grid::utils::idx3d(x, y, z, m_dims);
    }

    template <typename T> cl::sycl::buffer<T, 3> toSyclGrid(const std::vector<T>& vecGrid) const {
        return cl::sycl::buffer<T, 3>(vecGrid.data(), cl::sycl::range<3>(m_sizeX, m_sizeY, m_sizeZ));
    }

    template <typename T> std::vector<T> fromSyclGrid(cl::sycl::buffer<T, 3>& syclGrid) const {
        auto vecGrid = std::vector<T>(m_sizeX * m_sizeY * m_sizeZ);
        const auto gridAccessor = syclGrid.template get_access<cl::sycl::access::mode::read>();
        for (std::size_t x = 0; x < m_sizeX; ++x) {
            for (std::size_t y = 0; y < m_sizeY; ++y) {
                for (std::size_t z = 0; z < m_sizeZ; ++z) {
                    vecGrid[idx3d(x, y, z)] = gridAccessor[cl::sycl::id<3>(x, y, z)];
                }
            }
        }
        return vecGrid;
    }
};

template <class ProblemType, class Fields, unsigned padding>
RungeKuttaSyclSolver<ProblemType, Fields, padding>::RungeKuttaSyclSolver(const ProblemType& problem,
                                                                         const unsigned rungeKuttaSteps)
    : m_queue([](const cl::sycl::exception_list exceptions) {
          try {
              for (const auto& e : exceptions) {
                  std::rethrow_exception(e);
              }
          } catch (const cl::sycl::exception& e) {
              std::cout << "Exception during reduction: " << e.what() << std::endl;
          }
      }),
      m_sizeX(problem.numberCells[Direction::DirX] + 2 * padding),
      m_sizeY(problem.numberCells[Direction::DirY] + 2 * padding),
      m_sizeZ(problem.numberCells[Direction::DirZ] + 2 * padding),
      m_dims(grid::utils::dimensions{ m_sizeX, m_sizeY, m_sizeZ }), m_grid(m_sizeX * m_sizeY * m_sizeZ),
      changeBuffer(m_sizeX * m_sizeY * m_sizeZ), gridSubstepBuffer(m_sizeX * m_sizeY * m_sizeZ),
      m_cflBuffer(m_sizeX * m_sizeY * m_sizeZ), rungeKuttaSteps(rungeKuttaSteps), problem(problem),
      timeDelta(problem.timeDelta), timeCurrent(problem.timeStart), timeEnd(problem.timeEnd) {}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::initialise() {
    auto syclGrid = toSyclGrid(m_grid);
    problem.initialiseGridSycl(m_queue, syclGrid);
    m_queue.wait_and_throw();
    BoundarySycl::applyAll(m_queue, syclGrid, problem);
    m_queue.wait_and_throw();
    m_grid = fromSyclGrid(syclGrid);
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::step() {
    m_cfl = 0.0;
    for (unsigned substep = 0; substep < rungeKuttaSteps; substep++) {
        prepareSubstep();
        computeSubstep();
        m_cfl = reduceCFL(m_cfl, m_cflBuffer);
        finaliseSubstep(substep);
    }
    timeCurrent += timeDelta;
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::prepareSubstep() {
    // Number of negative temperatures is printed here, seems unnecessary

    // Start clock(s) -- time measurement omitted for now

    for (unsigned x = padding; x < m_sizeX - padding; ++x) {
        for (unsigned y = padding; y < m_sizeY - padding; ++y) {
            for (unsigned z = padding; z < m_sizeZ - padding; ++z) {
                using buffer_type = std::remove_cv_t<std::remove_reference_t<decltype(changeBuffer)>>;
                using value_type = typename buffer_type::value_type;
                changeBuffer[idx3d(x, y, z)] = value_type();
            }
        }
    }

    // CarbuncleFlag computation (not included by default)
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::computeSubstep() {
    for (unsigned x = padding; x < m_sizeX - padding + 1; ++x) {
        for (unsigned y = padding; y < m_sizeY - padding + 1; ++y) {
            for (unsigned z = padding; z < m_sizeZ - padding + 1; ++z) {
                const auto idx = idx3d(x, y, z);
                const auto item = cl::sycl::id<3>(x, y, z);
                std::tie(changeBuffer[idx], m_cflBuffer[idx]) = computeChanges(item);
            }
        }
    }
}

template <class ProblemType, class Fields, unsigned padding>
std::pair<Changes<Fields>, double>
RungeKuttaSyclSolver<ProblemType, Fields, padding>::computeChanges(const cl::sycl::id<3>& id) const {
    PerFaceValues reconstruction = ReconstructionSycl::reconstruct(m_grid, m_dims, id[0], id[1], id[2]);

    std::array<PhysValues, Faces::FaceMax> physicalValues;
    for (unsigned face = 0; face < Faces::FaceMax; face++) {
        Transformation::reconstToConservatives(physicalValues[face], reconstruction[face], problem.thermal,
                                               problem.gamma);
        physicalValues[face].thermalPressure =
            Transformation::computeThermalPressure(reconstruction[face], problem.thermal, problem.gamma);
        RiemannSolver::computeFluxes(physicalValues[face], reconstruction[face], face);
    }

    Changes<Fields> changes;
    auto localCFL = 0.0;
    for (unsigned dir = 0; dir < Direction::DirMax; dir++) {
        unsigned face = dir * 2;
        std::pair<double, double> characVelocities =
            RiemannSolver::characteristicVelocity(physicalValues[face], physicalValues[face + 1], reconstruction[face],
                                                  reconstruction[face + 1], problem.gamma, dir);
        localCFL = std::max(localCFL, calcCFL(characVelocities, dir));
        changes[dir] = RiemannSolver::numericalFlux(characVelocities, physicalValues[face], physicalValues[face + 1],
                                                    reconstruction[face], reconstruction[face + 1], dir);
    }
    return { changes, localCFL };
}

template <class ProblemType, class Fields, unsigned padding>
double RungeKuttaSyclSolver<ProblemType, Fields, padding>::calcCFL(const std::pair<double, double> characVelocities,
                                                                   const unsigned direction) const {
    const auto maxVelocity = std::max(characVelocities.first, characVelocities.second);
    const auto localCFL = maxVelocity * problem.inverseCellSize[direction];
    return localCFL;
}

template <class ProblemType, class Fields, unsigned padding>
double RungeKuttaSyclSolver<ProblemType, Fields, padding>::reduceCFL(const double initialCFL,
                                                                     const std::vector<double>& cflBuffer) const {
    auto cfl = initialCFL;
    for (unsigned x = padding; x < m_sizeX - padding + 1; ++x) {
        for (unsigned y = padding; y < m_sizeY - padding + 1; ++y) {
            for (unsigned z = padding; z < m_sizeZ - padding + 1; ++z) {
                const auto idx = idx3d(x, y, z);
                cfl = std::max(cfl, cflBuffer[idx]);
            }
        }
    }
    return cfl;
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::finaliseSubstep(const unsigned substep) {
    checkErrors();

    auto syclGrid = toSyclGrid(m_grid);
    problem.applySourceSycl(m_queue, syclGrid);
    m_grid = fromSyclGrid(syclGrid);
    checkErrors();

    syclGrid = toSyclGrid(m_grid);
    TransformationSycl::primitiveToConservative(m_queue, syclGrid, problem.thermal, problem.gamma);
    m_queue.wait_and_throw();
    m_grid = fromSyclGrid(syclGrid);
    integrateTime(substep);
    syclGrid = toSyclGrid(m_grid);
    TransformationSycl::conservativeToPrimitive(m_queue, syclGrid, problem.thermal, problem.gamma);
    m_queue.wait_and_throw();

    BoundarySycl::applyAll(m_queue, syclGrid, problem);
    m_queue.wait_and_throw();
    m_grid = fromSyclGrid(syclGrid);

    // Stop clock(s)
    // runtime estimation -- time measurement omitted for now

    // phystest
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::integrateTime(const unsigned substep) {
    if (substep == 0 && rungeKuttaSteps > 1) {
        gridSubstepBuffer = m_grid;
    }
    for (unsigned x = padding; x < m_sizeX - padding; x++) {
        for (unsigned y = padding; y < m_sizeY - padding; y++) {
            for (unsigned z = padding; z < m_sizeZ - padding; z++) {
                for (unsigned field = 0; field < Fields().size(); field++) {
                    const double changeX = (changeBuffer[idx3d(x + 1, y, z)][Direction::DirX][field] -
                                            changeBuffer[idx3d(x, y, z)][Direction::DirX][field]) *
                                           problem.inverseCellSize[Direction::DirX];
                    const double changeY = (changeBuffer[idx3d(x, y + 1, z)][Direction::DirY][field] -
                                            changeBuffer[idx3d(x, y, z)][Direction::DirY][field]) *
                                           problem.inverseCellSize[Direction::DirY];
                    const double changeZ = (changeBuffer[idx3d(x, y, z + 1)][Direction::DirZ][field] -
                                            changeBuffer[idx3d(x, y, z)][Direction::DirZ][field]) *
                                           problem.inverseCellSize[Direction::DirZ];
                    const double change = changeX + changeY + changeZ;
                    if (substep == 0) { // if-statement should be pulled out by compiler
                        m_grid[idx3d(x, y, z)][field] -= timeDelta * change;
                    } else {
                        m_grid[idx3d(x, y, z)][field] = 0.5 * gridSubstepBuffer[idx3d(x, y, z)][field] +
                                                        0.5 * m_grid[idx3d(x, y, z)][field] - 0.5 * timeDelta * change;
                    }
                }
            }
        }
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::adjust() {
    timeDelta = problem.cflThreshold / m_cfl;
    if (problem.preciseEnd) {
        timeDelta = std::min(timeDelta, timeEnd - timeCurrent);
    }
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::checkErrors() {
    if (GridFunctionsSycl::checkNaN(m_grid, m_dims)) {
        std::cerr << "Encountered NaN" << std::endl;
    }
}
