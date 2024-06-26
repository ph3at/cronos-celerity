#define CATCH_CONFIG_MAIN
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <toml++/toml.h>

#include <celerity.h>
#include <sycl/sycl.hpp>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

#include "../src/configuration/shock-tube.h"
#include "../src/grid/grid-functions.h"
#include "../src/solver/runge-kutta-celerity-solver.h"
#include "../src/solver/runge-kutta-solver.h"
#include "../src/solver/runge-kutta-sycl-solver.h"
#include "../src/sycl/reduction.h"

// This fixture (or a subclass) must be used by all tests that transitively use MPI.
class mpi_fixture {
  public:
    mpi_fixture() { celerity::detail::runtime::test_require_mpi(); }

    mpi_fixture(const mpi_fixture&) = delete;
    mpi_fixture& operator=(const mpi_fixture&) = delete;
};

// This fixture (or a subclass) must be used by all tests that transitively instantiate the runtime.
class runtime_fixture : public mpi_fixture {
  public:
    runtime_fixture() { celerity::detail::runtime::test_case_enter(); }

    runtime_fixture(const runtime_fixture&) = delete;
    runtime_fixture& operator=(const runtime_fixture&) = delete;

    ~runtime_fixture() {
        if (!celerity::detail::runtime::test_runtime_was_instantiated()) {
            WARN("Test specified a runtime_fixture, but did not end up instantiating the runtime");
        }
        celerity::detail::runtime::test_case_exit();
    }
};

class ScopedTimer {
  public:
    ScopedTimer(const std::string_view name) : m_name(name), m_start(std::chrono::high_resolution_clock::now()) {}
    ~ScopedTimer() {
        const auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - m_start).count() / 1000.0;
        std::cout << m_name << " took " << elapsed << "ms" << std::endl;
    }

  private:
    const std::string m_name;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};

TEST_CASE("Shock-Tube integration test", "[IntegrationTest]") {
    const toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");
    ShockTube shockTube(config);

    RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS> solver(shockTube);
    solver.initialise();

    double deviationPerStep = 1.0005;
    double deviationThreshold = deviationPerStep;
    for (unsigned timeStep = 1; timeStep <= 16; timeStep++) {
        solver.step();
        solver.adjust();
        const SimpleGrid<FieldStruct> baseline =
            GridFunctions::readFromFile("test-data/step-" + std::to_string(timeStep) + ".dat");
        double averageDeviation = GridFunctions::compare(baseline, solver.grid, false, false);
        REQUIRE(averageDeviation < deviationThreshold - 1.0);
        deviationThreshold *= deviationPerStep;
    }
}

TEST_CASE("Shock-Tube integration test with sycl", "[IntegrationTest][sycl]") {
    const toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");
    ShockTube shockTube(config);

    RungeKuttaSyclSolver<ShockTube, FieldStruct, GHOST_CELLS> solver(shockTube);
    solver.initialise();

    double deviationPerStep = 1.0005;
    double deviationThreshold = deviationPerStep;
    for (unsigned timeStep = 1; timeStep <= 16; timeStep++) {
        solver.step();
        solver.adjust();
        const SimpleGrid<FieldStruct> baseline =
            GridFunctions::readFromFile("test-data/step-" + std::to_string(timeStep) + ".dat");
        double averageDeviation = GridFunctions::compare(baseline, solver.grid(), false, false);
        REQUIRE(averageDeviation < deviationThreshold - 1.0);
        deviationThreshold *= deviationPerStep;
    }
}

TEST_CASE_METHOD(runtime_fixture, "Shock-Tube integration test with celerity", "[IntegrationTest][celerity]") {
    const toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");
    ShockTube shockTube(config);

    auto queue = celerity::distr_queue();
    RungeKuttaCeleritySolver<ShockTube, FieldStruct, GHOST_CELLS> solver(shockTube, queue);
    solver.initialise();

    double deviationPerStep = 1.0005;
    double deviationThreshold = deviationPerStep;
    for (unsigned timeStep = 1; timeStep <= 16; timeStep++) {
        solver.step();
        solver.adjust();
        const SimpleGrid<FieldStruct> baseline =
            GridFunctions::readFromFile("test-data/step-" + std::to_string(timeStep) + ".dat");
        double averageDeviation = GridFunctions::compare(baseline, solver.grid(), false, false);
        REQUIRE(averageDeviation < deviationThreshold - 1.0);
        deviationThreshold *= deviationPerStep;
    }
}

TEST_CASE("Shock-Tube integration test comparison host v sycl", "[IntegrationTest][sycl]") {
    const toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");
    const ShockTube shockTube(config);

    RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS> solver(shockTube);
    RungeKuttaSyclSolver<ShockTube, FieldStruct, GHOST_CELLS> syclSolver(shockTube);
    solver.initialise();
    syclSolver.initialise();

    double deviationPerStep = 1.0005;
    double deviationThreshold = deviationPerStep;
    while (!solver.isFinished() && !syclSolver.isFinished()) {
        solver.step();
        syclSolver.step();
        solver.adjust();
        syclSolver.adjust();
        const auto& baseline = solver.grid;
        double averageDeviation = GridFunctions::compare(baseline, syclSolver.grid(), false, false);
        // CHECK(averageDeviation == 0);
        REQUIRE(averageDeviation < deviationThreshold - 1.0);
        deviationThreshold *= deviationPerStep;
    }
    CHECK(solver.isFinished());
    CHECK(syclSolver.isFinished());
}

TEST_CASE_METHOD(runtime_fixture, "Shock-Tube integration test comparison host v celerity",
                 "[IntegrationTest][celerity]") {
    const toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");
    const ShockTube shockTube(config);

    RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS> solver(shockTube);
    auto queue = celerity::distr_queue();
    RungeKuttaCeleritySolver<ShockTube, FieldStruct, GHOST_CELLS> celeritySolver(shockTube, queue);
    solver.initialise();
    celeritySolver.initialise();

    double deviationPerStep = 1.0005;
    double deviationThreshold = deviationPerStep;
    while (!solver.isFinished() && !celeritySolver.isFinished()) {
        solver.step();
        celeritySolver.step();
        solver.adjust();
        celeritySolver.adjust();
        const auto& baseline = solver.grid;
        double averageDeviation = GridFunctions::compare(baseline, celeritySolver.grid(), false, false);
        // CHECK(averageDeviation == 0);
        REQUIRE(averageDeviation < deviationThreshold - 1.0);
        deviationThreshold *= deviationPerStep;
    }
    CHECK(solver.isFinished());
    CHECK(celeritySolver.isFinished());
}

TEST_CASE("Index mapping", "[sycl]") {
    SECTION("1D") {
        constexpr auto NUM_ELEMENTS = 1000;
        auto range = sycl::range<1>(NUM_ELEMENTS);
        auto data = std::vector<std::size_t>(range.size());
        for (std::size_t i = 0; i < data.size(); ++i) {
            data[i] = i;
        }
        auto buffer = sycl::buffer<std::size_t, 1>(data.data(), range);
        auto hostAccessor = buffer.get_host_access();
        for (std::size_t i = 0; i < data.size(); ++i) {
            const auto idx = sycl::id<1>(i);
            const auto dataElem = data[i];
            const auto bufferElem = hostAccessor[idx];
            CHECK(dataElem == bufferElem);
        }
    }
    SECTION("2D") {
        constexpr auto NUM_ELEMENTS = std::array<std::size_t, 2>{ 100, 300 };
        auto range = sycl::range<2>(NUM_ELEMENTS[0], NUM_ELEMENTS[1]);
        auto data = std::vector<std::size_t>(range.size());
        for (std::size_t i = 0; i < data.size(); ++i) {
            data[i] = i;
        }
        auto buffer = sycl::buffer<std::size_t, 2>(data.data(), range);
        auto hostAccessor = buffer.get_host_access();
        for (std::size_t i = 0; i < data.size(); ++i) {
            const auto idx1 = i % range[1];
            const auto idx0 = i / range[1];
            const auto idx = sycl::id<2>(idx0, idx1);
            const auto dataElem = data[i];
            const auto bufferElem = hostAccessor[idx];
            CHECK(dataElem == bufferElem);
        }
    }
    SECTION("3D") {
        constexpr auto NUM_ELEMENTS = std::array<std::size_t, 3>{ 500, 30, 100 };
        auto range = sycl::range<3>(NUM_ELEMENTS[0], NUM_ELEMENTS[1], NUM_ELEMENTS[2]);
        auto data = std::vector<std::size_t>(range.size());
        for (std::size_t i = 0; i < data.size(); ++i) {
            data[i] = i;
        }
        auto buffer = sycl::buffer<std::size_t, 3>(data.data(), range);
        auto hostAccessor = buffer.get_host_access();
        for (std::size_t i = 0; i < data.size(); ++i) {
            const auto idx2 = i % range[2];
            const auto idx1 = (i / range[2]) % range[1];
            const auto idx0 = i / range[2] / range[1];
            const auto idx = sycl::id<3>(idx0, idx1, idx2);
            const auto dataElem = data[i];
            const auto bufferElem = hostAccessor[idx];
            CHECK(dataElem == bufferElem);
        }
    }
}

TEST_CASE("Sycl", "[sycl]") {
    constexpr auto GHOST_CELLS = 2;
    constexpr int NUM_ELEMENTS[] = { 2, 3, 4 };
    constexpr int DATA_SIZE[] = {
        NUM_ELEMENTS[0] + 2 * GHOST_CELLS,
        NUM_ELEMENTS[1] + 2 * GHOST_CELLS,
        NUM_ELEMENTS[2] + 2 * GHOST_CELLS,
    };
    constexpr auto TOTAL_SIZE = DATA_SIZE[0] * DATA_SIZE[1] * DATA_SIZE[2];

    [[maybe_unused]] constexpr auto printData = [](const auto& data, const bool printGhostCells = true) {
        const auto ghost_cells = printGhostCells ? 0 : GHOST_CELLS;
        for (int d0 = ghost_cells; d0 < DATA_SIZE[0] - ghost_cells; ++d0) {
            for (int d1 = ghost_cells; d1 < DATA_SIZE[1] - ghost_cells; ++d1) {
                for (int d2 = ghost_cells; d2 < DATA_SIZE[2] - ghost_cells; ++d2) {
                    const auto idx = d0 * DATA_SIZE[1] * DATA_SIZE[2] + d1 * DATA_SIZE[2] + d2;
                    std::cout << data[idx] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    };

    auto dataA = std::vector<int>(TOTAL_SIZE);
    auto dataB = std::vector<int>(TOTAL_SIZE);
    auto dataC = std::vector<int>(TOTAL_SIZE, -1);

    for (std::size_t i = 0; i < TOTAL_SIZE; ++i) {
        const auto value = i;
        dataA[i] = value + 1;
        dataB[i] = 2 * (value + 1);
    }

    for (const auto& elem : dataC) {
        CHECK(elem == -1);
    }

    {
        auto queue = sycl::queue();
        auto itemRange = sycl::range<3>(NUM_ELEMENTS[0], NUM_ELEMENTS[1], NUM_ELEMENTS[2]);
        auto totalRange = sycl::range<3>(DATA_SIZE[0], DATA_SIZE[1], DATA_SIZE[2]);
        auto bufferA = sycl::buffer<int, 3>(dataA.data(), totalRange);
        auto bufferB = sycl::buffer<int, 3>(dataB.data(), totalRange);
        auto bufferC = sycl::buffer<int, 3>(dataC.data(), totalRange);

        queue.submit([&](sycl::handler& cgh) {
            auto accA = bufferA.template get_access<sycl::access::mode::read>(cgh);
            auto accB = bufferB.template get_access<sycl::access::mode::read>(cgh);
            auto accC = bufferC.template get_access<sycl::access::mode::discard_write>(cgh);

            cgh.parallel_for(itemRange, [=](const sycl::id<3> id) {
                const auto offset = sycl::id<3>(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);
                const auto i = id + offset;

                auto& res = accC[i];

                for (int d = 0; d < 3; ++d) {
                    for (int o = -GHOST_CELLS; o <= GHOST_CELLS; ++o) {
                        auto idx = i;
                        idx[d] += o;
                        res += accA[idx] + accB[idx];
                    }
                }
            });
        });
    }

    /*std::cout << "dataA:" << std::endl;
    printData(dataA);

    std::cout << "dataB:" << std::endl;
    printData(dataB);

    std::cout << "dataC:" << std::endl;
    printData(dataC);*/

    for (std::size_t d0 = GHOST_CELLS; d0 < DATA_SIZE[0] - GHOST_CELLS; ++d0) {
        for (std::size_t d1 = GHOST_CELLS; d1 < DATA_SIZE[1] - GHOST_CELLS; ++d1) {
            for (std::size_t d2 = GHOST_CELLS; d2 < DATA_SIZE[2] - GHOST_CELLS; ++d2) {
                const auto idx = d0 * DATA_SIZE[1] * DATA_SIZE[2] + d1 * DATA_SIZE[2] + d2;
                CHECK(dataC[idx] > 0);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

template <typename T> sycl::buffer<T, 1> uploadBuffer(const std::vector<T>& data) {
    return sycl::buffer<T, 1>(data.data(), sycl::range<1>(data.size()));
}

TEST_CASE("Sycl reduction", "[sycl]") {
    constexpr auto timedCreateQueue = []() {
        const auto _ = ScopedTimer("Creating SYCL queue");
        return sycl::queue([](const sycl::exception_list exceptions) {
            try {
                for (const auto& e : exceptions) {
                    std::rethrow_exception(e);
                }
            } catch (const sycl::exception& e) {
                std::cout << "Exception during reduction: " << e.what() << std::endl;
            }
        });
    };

    auto queue = timedCreateQueue();

    const auto NUM_ELEMENTS = GENERATE(1, 2, 3, 4, 128, 2048, 4096, 7757 /*, 512 * 1024 * 1024*/);
    auto data = std::vector<int>(NUM_ELEMENTS);

    {
        const auto _ = ScopedTimer("Generating " + std::to_string(NUM_ELEMENTS) + " random elements");
        std::iota(data.begin(), data.end(), -static_cast<int>(data.size()));
    }

    auto syclResult = 0;
    auto stlResult = 0;

    constexpr auto reductionOp = [](const int a, const int b) { return std::max(a, b); };

    {
        constexpr auto timedBufferUpload = [](const auto& data) {
            const auto _ = ScopedTimer("Uploading buffer to GPU");
            return uploadBuffer(data);
        };

        auto buffer = timedBufferUpload(data);

        {
            const auto _ = ScopedTimer("Sycl reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
            syclResult = sycl_utils::reduce(queue, buffer, reductionOp, std::numeric_limits<int>::lowest());
        }
    }
    {
        const auto _ = ScopedTimer("STL reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
        stlResult = std::accumulate(data.begin(), data.end(), std::numeric_limits<int>::lowest(), reductionOp);
    }

    CHECK(syclResult == stlResult);
    std::cout << std::endl;
}

#ifndef __HIPSYCL__

TEST_CASE("Sycl native reduction", "[sycl_reduction]") {
    constexpr auto timedCreateQueue = []() {
        const auto _ = ScopedTimer("Creating SYCL queue");
        return sycl::queue([](const sycl::exception_list exceptions) {
            try {
                for (const auto& e : exceptions) {
                    std::rethrow_exception(e);
                }
            } catch (const sycl::exception& e) {
                std::cout << "Exception during reduction: " << e.what() << std::endl;
            }
        });
    };

    auto queue = timedCreateQueue();

    const auto NUM_ELEMENTS = GENERATE(1, 2, 3, 4, 128, 2048, 4096, 7757 /*, 512 * 1024 * 1024*/);
    CAPTURE(NUM_ELEMENTS);
    auto data = std::vector<int>(NUM_ELEMENTS);

    {
        const auto _ = ScopedTimer("Generating " + std::to_string(NUM_ELEMENTS) + " random elements");
        std::iota(data.begin(), data.end(), -static_cast<int>(data.size()));
    }

    auto syclResult = 0;
    auto stlResult = 0;

    constexpr auto reductionOp = [](const int a, const int b) { return std::max(a, b); };
    {
        const auto _ = ScopedTimer("STL reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
        stlResult = std::accumulate(data.begin(), data.end(), std::numeric_limits<int>::lowest(), reductionOp);
    }
    {
        constexpr auto timedBufferUpload = [](const auto& data) {
            const auto _ = ScopedTimer("Uploading buffer to GPU");
            return uploadBuffer(data);
        };

        auto buffer = timedBufferUpload(data);

        auto maxBuffer = sycl::buffer<int, 1>{ { 1 } };

        {
            const auto _ = ScopedTimer("Celerity reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
            queue.submit([&](sycl::handler& cgh) {
                auto bufferAccessor = sycl::accessor{ buffer, cgh, sycl::read_only };

                auto maxReduction = sycl::reduction(maxBuffer, cgh, sycl::maximum<int>(),
                                                    sycl::property::reduction::initialize_to_identity{});

                cgh.parallel_for(buffer.get_range(), maxReduction,
                                 [=](sycl::item<1> idx, auto& max) { max.combine(bufferAccessor[idx]); });
            });
        }

        auto resultAccessor = sycl::host_accessor{ maxBuffer, sycl::read_only };
        syclResult = resultAccessor[0];
    }

    CHECK(syclResult == stlResult);
}

#endif

template <typename T> celerity::buffer<T, 1> uploadCelerityBuffer(const std::vector<T>& data) {
    return celerity::buffer<T, 1>(data.data(), celerity::range<1>(data.size()));
}

TEST_CASE_METHOD(runtime_fixture, "Celerity reduction", "[celerity]") {
    constexpr auto timedCreateQueue = []() {
        const auto _ = ScopedTimer("Creating celerity queue");
        return celerity::distr_queue();
    };

    auto queue = timedCreateQueue();

    const auto NUM_ELEMENTS = GENERATE(1, 2, 3, 4, 128, 2048, 4096, 7757 /*, 512 * 1024 * 1024*/);
    auto data = std::vector<int>(NUM_ELEMENTS);

    {
        const auto _ = ScopedTimer("Generating " + std::to_string(NUM_ELEMENTS) + " random elements");
        std::iota(data.begin(), data.end(), -static_cast<int>(data.size()));
    }

    auto celerityResult = 0;
    auto stlResult = 0;

    constexpr auto reductionOp = [](const int a, const int b) { return std::max(a, b); };
    {
        const auto _ = ScopedTimer("STL reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
        stlResult = std::accumulate(data.begin(), data.end(), std::numeric_limits<int>::lowest(), reductionOp);
    }
    {
        constexpr auto timedBufferUpload = [](const auto& data) {
            const auto _ = ScopedTimer("Uploading buffer to GPU");
            return uploadCelerityBuffer(data);
        };

        auto buffer = timedBufferUpload(data);

        auto maxBuffer = celerity::buffer<int, 1>{ { 1 } };

        {
            const auto _ = ScopedTimer("Celerity reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
            queue.submit([&buffer, &maxBuffer](celerity::handler& cgh) {
                auto bufferAccessor =
                    celerity::accessor{ buffer, cgh, celerity::access::one_to_one{}, celerity::read_only };

                auto maxReduction = celerity::reduction(maxBuffer, cgh, sycl::maximum<int>(),
                                                        celerity::property::reduction::initialize_to_identity{});

                celerity::debug::set_task_name(cgh, "reduction");

                cgh.parallel_for(buffer.get_range(), maxReduction,
                                 [=](celerity::item<1> idx, auto& max) { max.combine(bufferAccessor[idx]); });
            });
        }

        queue.submit([&maxBuffer, &celerityResult, &stlResult](celerity::handler& cgh) {
            celerity::accessor bufferAccessor{ maxBuffer, cgh, celerity::access::all{}, celerity::read_only_host_task };
            cgh.host_task(celerity::on_master_node, [bufferAccessor, &celerityResult, &stlResult] {
                celerityResult = bufferAccessor[0];
                CHECK(celerityResult == stlResult);
            });
        });

        queue.slow_full_sync();
    }
}

PaddedGrid<FieldStruct, GHOST_CELLS> syclGridToPaddedGrid(sycl::buffer<FieldStruct, 3>& grid) {
    const auto sizeX = grid.get_range()[0];
    const auto sizeY = grid.get_range()[1];
    const auto sizeZ = grid.get_range()[2];
    auto paddedGrid = PaddedGrid<FieldStruct, GHOST_CELLS>({}, sizeX - 2 * GHOST_CELLS, sizeY - 2 * GHOST_CELLS,
                                                           sizeZ - 2 * GHOST_CELLS);
    auto gridAccessor = sycl::host_accessor{ grid };
    for (std::size_t x = 0; x < sizeX; ++x) {
        for (std::size_t y = 0; y < sizeY; ++y) {
            for (std::size_t z = 0; z < sizeZ; ++z) {
                paddedGrid(x, y, z) = gridAccessor[sycl::id<3>(x, y, z)];
            }
        }
    }
    return paddedGrid;
}

PaddedGrid<FieldStruct, GHOST_CELLS> celerityGridToPaddedGrid(celerity::distr_queue& queue,
                                                              celerity::buffer<FieldStruct, 3>& grid) {
    const auto sizeX = grid.get_range()[0];
    const auto sizeY = grid.get_range()[1];
    const auto sizeZ = grid.get_range()[2];

    auto paddedGrid = PaddedGrid<FieldStruct, GHOST_CELLS>({}, sizeX - 2 * GHOST_CELLS, sizeY - 2 * GHOST_CELLS,
                                                           sizeZ - 2 * GHOST_CELLS);
    queue.submit([&paddedGrid, &grid, &sizeX, &sizeY, &sizeZ](celerity::handler& cgh) {
        auto gridAccessor = celerity::accessor{ grid, cgh, celerity::access::all{}, celerity::read_only_host_task };
        cgh.host_task(celerity::experimental::collective,
                      [=, &paddedGrid](celerity::experimental::collective_partition) {
                          for (std::size_t x = 0; x < sizeX; ++x) {
                              for (std::size_t y = 0; y < sizeY; ++y) {
                                  for (std::size_t z = 0; z < sizeZ; ++z) {
                                      paddedGrid(x, y, z) = gridAccessor[celerity::id<3>(x, y, z)];
                                  }
                              }
                          }
                      });
    });
    queue.slow_full_sync();
    return paddedGrid;
}

void compareGrids(const PaddedGrid<FieldStruct, GHOST_CELLS>& first,
                  const PaddedGrid<FieldStruct, GHOST_CELLS>& second) {
    REQUIRE(first.xDim() == second.xDim());
    REQUIRE(first.yDim() == second.yDim());
    REQUIRE(first.zDim() == second.zDim());
    REQUIRE(first.xStart() == second.xStart());
    REQUIRE(first.yStart() == second.yStart());
    REQUIRE(first.zStart() == second.zStart());
    REQUIRE(first.xEnd() == second.xEnd());
    REQUIRE(first.yEnd() == second.yEnd());
    REQUIRE(first.zEnd() == second.zEnd());

    auto failed = false;

    for (std::size_t x = 0; x < first.xDim(); ++x) {
        for (std::size_t y = 0; y < first.yDim(); ++y) {
            for (std::size_t z = 0; z < first.zDim(); ++z) {
                std::cout << "(";
                for (auto field = 0u; field < NUM_PHYSICAL_FIELDS; ++field) {
                    // CAPTURE(x);
                    // CAPTURE(y);
                    // CAPTURE(z);
                    // CAPTURE(field);
                    // CHECK(first(x, y, z)[field] == second(x, y, z)[field]);
                    failed |= first(x, y, z)[field] != second(x, y, z)[field];
                    const auto firstVal = first(x, y, z)[field];
                    const auto secondVal = second(x, y, z)[field];
                    // const auto diff = firstVal - secondVal;
                    // std::printf("%3.0f ", diff);
                    const auto* const res = firstVal < secondVal ? "<" : (firstVal > secondVal ? ">" : " ");
                    std::cout << res;
                }
                std::cout << ")";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    CHECK_FALSE(failed);
}

TEST_CASE("Sycl boundary v host", "[sycl][boundary]") {
    auto queue = sycl::queue();

    toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");

    *(config["grid"]["number_cells_x"].as_integer()) = 5;
    *(config["grid"]["number_cells_y"].as_integer()) = 6;
    *(config["grid"]["number_cells_z"].as_integer()) = 7;

    const auto BOUNDARY = GENERATE("outflow", "user", "extrapolate");
    *(config["boundary"]["west"].as_string()) = BOUNDARY;
    *(config["boundary"]["east"].as_string()) = BOUNDARY;
    *(config["boundary"]["south"].as_string()) = BOUNDARY;
    *(config["boundary"]["north"].as_string()) = BOUNDARY;
    *(config["boundary"]["bottom"].as_string()) = BOUNDARY;
    *(config["boundary"]["top"].as_string()) = BOUNDARY;

    const ShockTube shockTube(config);

    const auto sizeX = shockTube.numberCells[Direction::DirX];
    const auto sizeY = shockTube.numberCells[Direction::DirY];
    const auto sizeZ = shockTube.numberCells[Direction::DirZ];

    auto grid = sycl::buffer<FieldStruct, 3>(sycl::range<3>(sizeX, sizeY, sizeZ));

    queue.submit([&](sycl::handler& cgh) {
        auto gridAcc = sycl::accessor{ grid, cgh, sycl::write_only, sycl::no_init };

        cgh.parallel_for(grid.get_range(), [=](const sycl::id<3> id) {
            const auto linId = id[0] * sizeY * sizeZ + id[1] * sizeZ + id[2] + 1;
            auto val = FieldStruct{};
            for (auto i = 0u; i < NUM_PHYSICAL_FIELDS; ++i) {
                val[i] = linId * NUM_PHYSICAL_FIELDS + i;
            }
            gridAcc[id] = val;
        });
    });

    auto paddedGrid = syclGridToPaddedGrid(grid);

    BoundarySycl::applyAll(queue, grid, shockTube);
    Boundary::applyAll(paddedGrid, shockTube);

    const auto syclBoundary = syclGridToPaddedGrid(grid);
    compareGrids(syclBoundary, paddedGrid);
}

TEST_CASE_METHOD(runtime_fixture, "Celerity boundary v host", "[celerity][boundary]") {
    auto queue = celerity::distr_queue();

    toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");

    *(config["grid"]["number_cells_x"].as_integer()) = 5;
    *(config["grid"]["number_cells_y"].as_integer()) = 6;
    *(config["grid"]["number_cells_z"].as_integer()) = 7;

    const auto BOUNDARY = GENERATE("outflow", "user", "extrapolate");
    *(config["boundary"]["west"].as_string()) = BOUNDARY;
    *(config["boundary"]["east"].as_string()) = BOUNDARY;
    *(config["boundary"]["south"].as_string()) = BOUNDARY;
    *(config["boundary"]["north"].as_string()) = BOUNDARY;
    *(config["boundary"]["bottom"].as_string()) = BOUNDARY;
    *(config["boundary"]["top"].as_string()) = BOUNDARY;

    const ShockTube shockTube(config);

    const auto sizeX = shockTube.numberCells[Direction::DirX];
    const auto sizeY = shockTube.numberCells[Direction::DirY];
    const auto sizeZ = shockTube.numberCells[Direction::DirZ];

    auto grid = celerity::buffer<FieldStruct, 3>(celerity::range<3>(sizeX, sizeY, sizeZ));

    queue.submit([&](celerity::handler& cgh) {
        auto gridAcc =
            celerity::accessor{ grid, cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init };

        celerity::debug::set_task_name(cgh, "testGridInit");

        cgh.parallel_for(grid.get_range(), [=](const celerity::id<3> id) {
            const auto linId = id[0] * sizeY * sizeZ + id[1] * sizeZ + id[2] + 1;
            auto val = FieldStruct{};
            for (auto i = 0u; i < NUM_PHYSICAL_FIELDS; ++i) {
                val[i] = linId * NUM_PHYSICAL_FIELDS + i;
            }
            gridAcc[id] = val;
        });
    });

    auto paddedGrid = celerityGridToPaddedGrid(queue, grid);

    BoundaryCelerity::applyAll(queue, grid, shockTube);
    Boundary::applyAll(paddedGrid, shockTube);

    const auto celerityBoundary = celerityGridToPaddedGrid(queue, grid);
    compareGrids(celerityBoundary, paddedGrid);
}

TEST_CASE_METHOD(runtime_fixture, "Celerity boundary v sycl", "[celerity][sycl][boundary]") {
    toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");

    *(config["grid"]["number_cells_x"].as_integer()) = 5;
    *(config["grid"]["number_cells_y"].as_integer()) = 6;
    *(config["grid"]["number_cells_z"].as_integer()) = 7;

    const auto BOUNDARY = GENERATE("outflow", "user", "extrapolate");
    *(config["boundary"]["west"].as_string()) = BOUNDARY;
    *(config["boundary"]["east"].as_string()) = BOUNDARY;
    *(config["boundary"]["south"].as_string()) = BOUNDARY;
    *(config["boundary"]["north"].as_string()) = BOUNDARY;
    *(config["boundary"]["bottom"].as_string()) = BOUNDARY;
    *(config["boundary"]["top"].as_string()) = BOUNDARY;

    const ShockTube shockTube(config);

    const auto sizeX = shockTube.numberCells[Direction::DirX];
    const auto sizeY = shockTube.numberCells[Direction::DirY];
    const auto sizeZ = shockTube.numberCells[Direction::DirZ];

    const auto celerityBoundary = [&]() {
        auto celerityGrid = celerity::buffer<FieldStruct, 3>(celerity::range<3>(sizeX, sizeY, sizeZ));
        auto distrQueue = celerity::distr_queue();
        distrQueue.submit([&](celerity::handler& cgh) {
            auto gridAcc = celerity::accessor{ celerityGrid, cgh, celerity::access::one_to_one{}, celerity::write_only,
                                               celerity::no_init };

            celerity::debug::set_task_name(cgh, "testGridInit");

            cgh.parallel_for(celerityGrid.get_range(), [=](const celerity::id<3> id) {
                const auto linId = id[0] * sizeY * sizeZ + id[1] * sizeZ + id[2] + 1;
                auto val = FieldStruct{};
                for (auto i = 0u; i < NUM_PHYSICAL_FIELDS; ++i) {
                    val[i] = linId * NUM_PHYSICAL_FIELDS + i;
                }
                gridAcc[id] = val;
            });
        });
        BoundaryCelerity::applyAll(distrQueue, celerityGrid, shockTube);

        return celerityGridToPaddedGrid(distrQueue, celerityGrid);
    }();

    const auto syclBoundary = [&]() {
        auto queue = sycl::queue();
        auto syclGrid = sycl::buffer<FieldStruct, 3>(sycl::range<3>(sizeX, sizeY, sizeZ));
        queue.submit([&](sycl::handler& cgh) {
            auto gridAcc = sycl::accessor{ syclGrid, cgh, sycl::write_only, sycl::no_init };

            cgh.parallel_for(syclGrid.get_range(), [=](const sycl::id<3> id) {
                const auto linId = id[0] * sizeY * sizeZ + id[1] * sizeZ + id[2] + 1;
                auto val = FieldStruct{};
                for (auto i = 0u; i < NUM_PHYSICAL_FIELDS; ++i) {
                    val[i] = linId * NUM_PHYSICAL_FIELDS + i;
                }
                gridAcc[id] = val;
            });
        });
        BoundarySycl::applyAll(queue, syclGrid, shockTube);
        return syclGridToPaddedGrid(syclGrid);
    }();

    compareGrids(celerityBoundary, syclBoundary);
}

TEST_CASE_METHOD(runtime_fixture, "3D celerity boundary v host", "[celerity][boundary]") {
    auto queue = celerity::distr_queue();

    toml::table config = toml::parse_file("configuration/shock-tube-integration.toml");

    *(config["grid"]["number_cells_x"].as_integer()) = 5;
    *(config["grid"]["number_cells_y"].as_integer()) = 6;
    *(config["grid"]["number_cells_z"].as_integer()) = 7;

    const auto BOUNDARY = GENERATE("outflow", "user", "extrapolate");
    *(config["boundary"]["west"].as_string()) = BOUNDARY;
    *(config["boundary"]["east"].as_string()) = BOUNDARY;
    *(config["boundary"]["south"].as_string()) = BOUNDARY;
    *(config["boundary"]["north"].as_string()) = BOUNDARY;
    *(config["boundary"]["bottom"].as_string()) = BOUNDARY;
    *(config["boundary"]["top"].as_string()) = BOUNDARY;

    const ShockTube shockTube(config);

    const auto sizeX = shockTube.numberCells[Direction::DirX];
    const auto sizeY = shockTube.numberCells[Direction::DirY];
    const auto sizeZ = shockTube.numberCells[Direction::DirZ];

    auto grid = celerity::buffer<FieldStruct, 3>(celerity::range<3>(sizeX, sizeY, sizeZ));

    queue.submit([&](celerity::handler& cgh) {
        auto gridAcc =
            celerity::accessor{ grid, cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init };

        celerity::debug::set_task_name(cgh, "testGridInit3D");

        cgh.parallel_for(grid.get_range(), [=](const celerity::id<3> id) {
            const auto linId = id[0] * sizeY * sizeZ + id[1] * sizeZ + id[2] + 1;
            auto val = FieldStruct{};
            for (auto i = 0u; i < NUM_PHYSICAL_FIELDS; ++i) {
                val[i] = linId * NUM_PHYSICAL_FIELDS + i;
            }
            gridAcc[id] = val;
        });
    });

    auto paddedGrid = celerityGridToPaddedGrid(queue, grid);

    BoundaryCelerity::applyAll3D(queue, grid, shockTube);
    Boundary::applyAll(paddedGrid, shockTube);

    const auto celerityBoundary = celerityGridToPaddedGrid(queue, grid);
    compareGrids(celerityBoundary, paddedGrid);
}

int main(int argc, char* argv[]) {
    Catch::Session session;

    int return_code = session.applyCommandLine(argc, argv);
    if (return_code != 0) {
        return return_code;
    }

    celerity::detail::runtime::test_mode_enter();
    return_code = session.run();
    celerity::detail::runtime::test_mode_exit();
    return return_code;
}
