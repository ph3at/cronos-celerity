#pragma once

#include <algorithm>
#include <array>
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
    void integrateTime(cl::sycl::queue& queue, cl::sycl::buffer<Fields, 3>& grid,
                       cl::sycl::buffer<Fields, 3>& gridSubstep, cl::sycl::buffer<Changes<Fields>, 3>& changes,
                       const unsigned substep) const;
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
    : m_queue(cl::sycl::gpu_selector{},
              [](const cl::sycl::exception_list exceptions) {
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

    TransformationSycl::primitiveToConservative(m_queue, syclGrid, problem.thermal, problem.gamma);
    m_queue.wait_and_throw();
    auto syclGridSubstep = toSyclGrid(gridSubstepBuffer);
    auto syclChanges = toSyclGrid(changeBuffer);
    integrateTime(m_queue, syclGrid, syclGridSubstep, syclChanges, substep);
    m_queue.wait_and_throw();
    gridSubstepBuffer = fromSyclGrid(syclGridSubstep);
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
void RungeKuttaSyclSolver<ProblemType, Fields, padding>::integrateTime(cl::sycl::queue& queue,
                                                                       cl::sycl::buffer<Fields, 3>& grid,
                                                                       cl::sycl::buffer<Fields, 3>& gridSubstep,
                                                                       cl::sycl::buffer<Changes<Fields>, 3>& changes,
                                                                       const unsigned substep) const {
    queue.submit([&](cl::sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<cl::sycl::access::mode::read_write>(cgh);
        auto gridSubstepAccessor = gridSubstep.template get_access<cl::sycl::access::mode::read_write>(cgh);
        auto changeAccessor = changes.template get_access<cl::sycl::access::mode::read>(cgh);

        auto range = grid.get_range();
        range[0] -= 2 * padding;
        range[1] -= 2 * padding;
        range[2] -= 2 * padding;

        cgh.parallel_for(range, [=, timeDelta = timeDelta, rungeKuttaSteps = rungeKuttaSteps,
                                 inverseCellSize = problem.inverseCellSize](const cl::sycl::id<3> id) {
            const auto offset = cl::sycl::id<3>(padding, padding, padding);
            const auto idx = id + offset;

            if (substep == 0 && rungeKuttaSteps > 1) {
                gridSubstepAccessor[idx] = gridAccessor[idx];
            }

            for (unsigned field = 0; field < std::tuple_size<Fields>{}; ++field) {
                const auto nextXIdx = cl::sycl::id<3>(idx[0] + 1, idx[1], idx[2]);
                const auto nextYIdx = cl::sycl::id<3>(idx[0], idx[1] + 1, idx[2]);
                const auto nextZIdx = cl::sycl::id<3>(idx[0], idx[1], idx[2] + 1);

                const double changeX =
                    (changeAccessor[nextXIdx][Direction::DirX][field] - changeAccessor[idx][Direction::DirX][field]) *
                    inverseCellSize[Direction::DirX];
                const double changeY =
                    (changeAccessor[nextYIdx][Direction::DirY][field] - changeAccessor[idx][Direction::DirY][field]) *
                    inverseCellSize[Direction::DirY];
                const double changeZ =
                    (changeAccessor[nextZIdx][Direction::DirZ][field] - changeAccessor[idx][Direction::DirZ][field]) *
                    inverseCellSize[Direction::DirZ];
                const double change = changeX + changeY + changeZ;
                if (substep == 0) {
                    gridAccessor[idx][field] -= timeDelta * change;
                } else {
                    gridAccessor[idx][field] = 0.5 * gridSubstepAccessor[idx][field] + 0.5 * gridAccessor[idx][field] -
                                               0.5 * timeDelta * change;
                }
            }
        });
    });
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
