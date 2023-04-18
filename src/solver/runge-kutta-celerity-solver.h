#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <tuple>
#include <vector>

#include <celerity.h>
#include <sycl/sycl.hpp>

#include "../boundary/boundary.h"
#include "../data-types/direction.h"
#include "../data-types/faces.h"
#include "../grid/grid-functions.h"
#include "../grid/padded-grid.h"
#include "../riemann/reconstruction.h"
#include "../riemann/riemann-solver.h"
#include "../sycl/reduction.h"
#include "../transformation/transformations.h"

template <class Fields> using Changes = std::array<Fields, Direction::DirMax>;

template <class ProblemType, class Fields, unsigned padding> class RungeKuttaCeleritySolver {
  public:
    RungeKuttaCeleritySolver(celerity::distr_queue& celerity_queue, const ProblemType& problem,
                             const unsigned rungeKuttaSteps = 2);

    void initialise();
    void step();
    void adjust();
    void finalise(){};
    bool isFinished() const;
    void report(const int timeSteps) const;

    PaddedGrid<Fields, padding> grid() {
        auto paddedGrid =
            PaddedGrid<Fields, padding>({}, m_sizeX - 2 * padding, m_sizeY - 2 * padding, m_sizeZ - 2 * padding);
        const auto gridAccessor = m_grid.template get_access<sycl::access::mode::read>();
        for (std::size_t x = 0; x < m_sizeX; ++x) {
            for (std::size_t y = 0; y < m_sizeY; ++y) {
                for (std::size_t z = 0; z < m_sizeZ; ++z) {
                    paddedGrid(x, y, z) = gridAccessor[sycl::id<3>(x, y, z)];
                }
            }
        }
        return paddedGrid;
    }

  private:
    sycl::queue m_queue;
    celerity::distr_queue& m_celerity_queue;

    const std::size_t m_sizeX;
    const std::size_t m_sizeY;
    const std::size_t m_sizeZ;

    sycl::buffer<Fields, 3> m_grid;
    sycl::buffer<Changes<Fields>, 3> m_changeBuffer;
    sycl::buffer<Fields, 3> m_gridSubstepBuffer;
    sycl::buffer<double, 1> m_cflBuffer;

    double m_cfl;
    const unsigned rungeKuttaSteps;

    const ProblemType& problem;

    double timeDelta;
    double timeCurrent;
    const double timeEnd;

    template <typename T, int Dims>
    celerity::buffer<T, Dims> convertSyclToCelerityBuffer(sycl::buffer<T, Dims>& buffer) const {
        auto data = std::vector<T>(buffer.size());
        auto hostAccessor = buffer.get_host_access();
        for (std::size_t i = 0; i < buffer.size(); ++i) {
            const auto id = idFromIdx(buffer.get_range(), i);
            data[i] = hostAccessor[id];
        }
        auto celerityBuffer = celerity::buffer<T, Dims>(data.data(), buffer.get_range());
        return celerityBuffer;
    }

    template <int Dims> celerity::id<Dims> idFromIdx(const celerity::range<Dims>& range, const std::size_t idx) const {
        if constexpr (Dims == 1) {
            return celerity::id<Dims>(idx);
        }
        if constexpr (Dims == 2) {
            const auto idx1 = idx % range[1];
            const auto idx0 = idx / range[1];
            return celerity::id<Dims>(idx0, idx1);
        }
        if constexpr (Dims == 3) {
            const auto idx2 = idx % range[2];
            const auto idx1 = (idx / range[2]) % range[1];
            const auto idx0 = idx / range[2] / range[1];
            return celerity::id<Dims>(idx0, idx1, idx2);
        }
    }

    template <typename T, int Dims>
    sycl::buffer<T, Dims> convertCelerityToSyclBuffer(celerity::buffer<T, Dims>& buffer) const {
        auto data = std::vector<T>(buffer.get_range().size());
        m_celerity_queue.submit(celerity::allow_by_ref, [=, &data](celerity::handler& cgh) {
            celerity::accessor bufferAccessor{ buffer, cgh, celerity::access::all{}, celerity::read_only_host_task };
            cgh.host_task(celerity::experimental::collective, [=, &data](celerity::experimental::collective_partition) {
                for (std::size_t i = 0; i < data.size(); ++i) {
                    const auto id = idFromIdx(buffer.get_range(), i);
                    data[i] = bufferAccessor[id];
                }
            });
        });
        m_celerity_queue.slow_full_sync();

        auto syclBuffer = sycl::buffer<T, Dims>(buffer.get_range());
        auto hostAccessor = syclBuffer.get_host_access();
        for (std::size_t i = 0; i < syclBuffer.size(); ++i) {
            const auto id = idFromIdx(buffer.get_range(), i);
            hostAccessor[id] = data[i];
        }
        return syclBuffer;
    }

    void prepareSubstep(celerity::distr_queue& queue, celerity::buffer<Changes<Fields>, 3>& changeBuffer) const;
    void computeSubstep();
    static double calcCFL(std::pair<double, double> characVelocities, const double inverseCellSize);
    double reduceCFL(celerity::distr_queue& queue, sycl::buffer<double, 1>& cflBuffer, const double initialCFL);
    void finaliseSubstep(const unsigned substep);
    void integrateTime(sycl::queue& queue, sycl::buffer<Fields, 3>& grid, sycl::buffer<Fields, 3>& gridSubstep,
                       sycl::buffer<Changes<Fields>, 3>& changes, const unsigned substep) const;
    void checkErrors(celerity::distr_queue& queue, celerity::buffer<Fields, 3>& grid) const;
};

template <class ProblemType, class Fields, unsigned padding>
RungeKuttaCeleritySolver<ProblemType, Fields, padding>::RungeKuttaCeleritySolver(celerity::distr_queue& celerity_queue,
                                                                                 const ProblemType& problem,
                                                                                 const unsigned rungeKuttaSteps)
    : m_queue(sycl::gpu_selector{},
              [](const sycl::exception_list exceptions) {
                  try {
                      for (const auto& e : exceptions) {
                          std::rethrow_exception(e);
                      }
                  } catch (const sycl::exception& e) {
                      std::cout << "Exception during reduction: " << e.what() << std::endl;
                  }
              }),
      m_celerity_queue(celerity_queue), m_sizeX(problem.numberCells[Direction::DirX] + 2 * padding),
      m_sizeY(problem.numberCells[Direction::DirY] + 2 * padding),
      m_sizeZ(problem.numberCells[Direction::DirZ] + 2 * padding), m_grid(sycl::range<3>(m_sizeX, m_sizeY, m_sizeZ)),
      m_changeBuffer(sycl::range<3>(m_sizeX, m_sizeY, m_sizeZ)),
      m_gridSubstepBuffer(sycl::range<3>(m_sizeX, m_sizeY, m_sizeZ)),
      m_cflBuffer(sycl::range<1>(m_sizeX * m_sizeY * m_sizeZ)), rungeKuttaSteps(rungeKuttaSteps), problem(problem),
      timeDelta(problem.timeDelta), timeCurrent(problem.timeStart), timeEnd(problem.timeEnd) {
    m_queue.submit([&](sycl::handler& cgh) {
        auto cflAccessor = m_cflBuffer.template get_access<sycl::access::mode::discard_write>(cgh);
        cgh.fill(cflAccessor, std::numeric_limits<double>::lowest());
    });
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::initialise() {
    auto celerityGrid = convertSyclToCelerityBuffer(m_grid);
    problem.initialiseGridCelerity(m_celerity_queue, celerityGrid);
    m_grid = convertCelerityToSyclBuffer(celerityGrid);
    BoundarySycl::applyAll(m_queue, m_grid, problem);
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::step() {
    m_cfl = 0.0;
    for (unsigned substep = 0; substep < rungeKuttaSteps; substep++) {
        auto celerityChangeBuffer = convertSyclToCelerityBuffer(m_changeBuffer);
        prepareSubstep(m_celerity_queue, celerityChangeBuffer);
        m_changeBuffer = convertCelerityToSyclBuffer(celerityChangeBuffer);
        computeSubstep();
        m_cfl = reduceCFL(m_celerity_queue, m_cflBuffer, m_cfl);
        finaliseSubstep(substep);
    }
    timeCurrent += timeDelta;
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::prepareSubstep(
    celerity::distr_queue& queue, celerity::buffer<Changes<Fields>, 3>& changeBuffer) const {
    // Number of negative temperatures is printed here, seems unnecessary

    // Start clock(s) -- time measurement omitted for now

    // TODO: Clearing the buffer is not necessary -> get rid of it

    queue.submit([=](celerity::handler& cgh) {
        auto changeAccessor = celerity::accessor{ changeBuffer, cgh, celerity::access::one_to_one{},
                                                  celerity::write_only, celerity::no_init };

        cgh.parallel_for(changeBuffer.get_range(),
                         [=](const celerity::id<3> id) { changeAccessor[id] = Changes<Fields>(); });
    });
    queue.slow_full_sync();

    // CarbuncleFlag computation (not included by default)
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::computeSubstep() {
    m_queue.submit([&](sycl::handler& cgh) {
        auto gridAccessor = m_grid.template get_access<sycl::access::mode::read>(cgh);
        auto changeAccessor = m_changeBuffer.template get_access<sycl::access::mode::write>(cgh);
        auto cflAccessor = m_cflBuffer.template get_access<sycl::access::mode::write>(cgh);

        auto range = m_grid.get_range();
        range[0] = range[0] - 2 * padding + 1;
        range[1] = range[1] - 2 * padding + 1;
        range[2] = range[2] - 2 * padding + 1;

        cgh.parallel_for(range, [=, thermal = problem.thermal, gamma = problem.gamma,
                                 inverseCellSize = problem.inverseCellSize](const sycl::id<3> id) {
            const auto offset = sycl::id<3>(padding, padding, padding);
            const auto idx = id + offset;

            PerFaceValues reconstruction = ReconstructionSycl::reconstruct(gridAccessor, idx);

            std::array<PhysValues, Faces::FaceMax> physicalValues;
            for (unsigned face = 0; face < Faces::FaceMax; face++) {
                Transformation::reconstToConservatives(physicalValues[face], reconstruction[face], thermal, gamma);
                physicalValues[face].thermalPressure =
                    Transformation::computeThermalPressure(reconstruction[face], thermal, gamma);
                RiemannSolver::computeFluxes(physicalValues[face], reconstruction[face], face);
            }

            Changes<Fields> changes;
            auto localCFL = 0.0;
            for (unsigned dir = 0; dir < Direction::DirMax; dir++) {
                unsigned face = dir * 2;
                std::pair<double, double> characVelocities =
                    RiemannSolver::characteristicVelocity(physicalValues[face], physicalValues[face + 1],
                                                          reconstruction[face], reconstruction[face + 1], gamma, dir);
                localCFL = std::max(localCFL, calcCFL(characVelocities, inverseCellSize[dir]));
                changes[dir] =
                    RiemannSolver::numericalFlux(characVelocities, physicalValues[face], physicalValues[face + 1],
                                                 reconstruction[face], reconstruction[face + 1], dir);
            }

            changeAccessor[idx] = changes;
            const auto linIdx = sycl::id<1>(idx[0] * range[1] * range[2] + idx[1] * range[2] + idx[2]);
            cflAccessor[linIdx] = localCFL;
        });
    });
}

template <class ProblemType, class Fields, unsigned padding>
double RungeKuttaCeleritySolver<ProblemType, Fields, padding>::calcCFL(const std::pair<double, double> characVelocities,
                                                                       const double inverseCellSize) {
    const auto maxVelocity = std::max(characVelocities.first, characVelocities.second);
    const auto localCFL = maxVelocity * inverseCellSize;
    return localCFL;
}

template <class ProblemType, class Fields, unsigned padding>
double RungeKuttaCeleritySolver<ProblemType, Fields, padding>::reduceCFL(celerity::distr_queue& queue,
                                                                         sycl::buffer<double, 1>& cflBuffer,
                                                                         const double initialCFL) {
    auto resultBuffer = celerity::buffer<double, 1>{ celerity::range<1>{ 1 } };
    auto celerityCflBuffer = convertSyclToCelerityBuffer(cflBuffer);

    queue.submit([=](celerity::handler& cgh) {
        auto bufferAccessor =
            celerity::accessor{ celerityCflBuffer, cgh, celerity::access::one_to_one{}, celerity::read_only };
        auto maxReduction = celerity::reduction(resultBuffer, cgh, sycl::maximum<>(),
                                                celerity::property::reduction::initialize_to_identity{});

        cgh.parallel_for(celerityCflBuffer.get_range(), maxReduction,
                         [=](celerity::item<1> idx, auto& max) { max.combine(bufferAccessor[idx]); });
    });

    double reducedCFL = 0;

    queue.submit(celerity::allow_by_ref, [=, &reducedCFL](celerity::handler& cgh) {
        celerity::accessor bufferAccessor{ resultBuffer, cgh, celerity::access::all{}, celerity::read_only_host_task };
        cgh.host_task(
            celerity::experimental::collective,
            [=, &reducedCFL](celerity::experimental::collective_partition) { reducedCFL = bufferAccessor[0]; });
    });

    queue.slow_full_sync();
    const auto maxCFL = std::max(initialCFL, reducedCFL);
    return maxCFL;
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::finaliseSubstep(const unsigned substep) {
#ifndef NDEBUG
    auto celerityGrid = convertSyclToCelerityBuffer(m_grid);
    checkErrors(m_celerity_queue, celerityGrid);
#endif

    problem.applySourceSycl(m_queue, m_grid);

#ifndef NDEBUG
    celerityGrid = convertSyclToCelerityBuffer(m_grid);
    checkErrors(m_celerity_queue, celerityGrid);
#endif

    TransformationSycl::primitiveToConservative(m_queue, m_grid, problem.thermal, problem.gamma);
    integrateTime(m_queue, m_grid, m_gridSubstepBuffer, m_changeBuffer, substep);
    TransformationSycl::conservativeToPrimitive(m_queue, m_grid, problem.thermal, problem.gamma);

    BoundarySycl::applyAll(m_queue, m_grid, problem);

    // Stop clock(s)
    // runtime estimation -- time measurement omitted for now

    // phystest
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::integrateTime(sycl::queue& queue,
                                                                           sycl::buffer<Fields, 3>& grid,
                                                                           sycl::buffer<Fields, 3>& gridSubstep,
                                                                           sycl::buffer<Changes<Fields>, 3>& changes,
                                                                           const unsigned substep) const {
    queue.submit([&](sycl::handler& cgh) {
        auto gridAccessor = grid.template get_access<sycl::access::mode::read_write>(cgh);
        auto gridSubstepAccessor = (substep == 0 && rungeKuttaSteps > 1)
                                       ? sycl::accessor(gridSubstep, cgh, sycl::read_write, sycl::no_init)
                                       : sycl::accessor(gridSubstep, cgh, sycl::read_write);
        auto changeAccessor = changes.template get_access<sycl::access::mode::read>(cgh);

        auto range = grid.get_range();
        range[0] -= 2 * padding;
        range[1] -= 2 * padding;
        range[2] -= 2 * padding;

        cgh.parallel_for(range, [=, timeDelta = timeDelta, rungeKuttaSteps = rungeKuttaSteps,
                                 inverseCellSize = problem.inverseCellSize](const sycl::id<3> id) {
            const auto offset = sycl::id<3>(padding, padding, padding);
            const auto idx = id + offset;

            if (substep == 0 && rungeKuttaSteps > 1) {
                gridSubstepAccessor[idx] = gridAccessor[idx];
            }

            for (unsigned field = 0; field < std::tuple_size<Fields>{}; ++field) {
                const auto nextXIdx = sycl::id<3>(idx[0] + 1, idx[1], idx[2]);
                const auto nextYIdx = sycl::id<3>(idx[0], idx[1] + 1, idx[2]);
                const auto nextZIdx = sycl::id<3>(idx[0], idx[1], idx[2] + 1);

                const auto calcDirChanges = [&](const auto& nextIdx, const auto& dir) {
                    return (changeAccessor[nextIdx][dir][field] - changeAccessor[idx][dir][field]) *
                           inverseCellSize[dir];
                };

                const auto changeX = calcDirChanges(nextXIdx, Direction::DirX);
                const auto changeY = calcDirChanges(nextYIdx, Direction::DirY);
                const auto changeZ = calcDirChanges(nextZIdx, Direction::DirZ);
                const auto change = changeX + changeY + changeZ;

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
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::adjust() {
    timeDelta = problem.cflThreshold / m_cfl;
    if (problem.preciseEnd) {
        timeDelta = std::min(timeDelta, timeEnd - timeCurrent);
    }
}

template <class ProblemType, class Fields, unsigned padding>
bool RungeKuttaCeleritySolver<ProblemType, Fields, padding>::isFinished() const {
    return timeCurrent >= timeEnd;
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::report(const int timeSteps) const {
    std::cout << "Stopped at time " << timeCurrent << ", after " << timeSteps << " steps." << std::endl;
}

template <class ProblemType, class Fields, unsigned padding>
void RungeKuttaCeleritySolver<ProblemType, Fields, padding>::checkErrors(celerity::distr_queue& queue,
                                                                         celerity::buffer<Fields, 3>& grid) const {
    if (GridFunctionsCelerity::checkNaN(queue, grid)) {
        std::cerr << "Encountered NaN" << std::endl;
    }
}
