#define CATCH_CONFIG_MAIN
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
#include "../src/solver/runge-kutta-solver.h"
#include "../src/solver/runge-kutta-sycl-solver.h"
#include "../src/sycl/reduction.h"

static celerity::distr_queue& get_celerity_queue() {
    static celerity::distr_queue queue;
    return queue;
}

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

template <typename Solver>
auto benchmarkSolver(const int numSizes, const int numRuns, const std::string_view configPath) {
    auto results = std::vector(numSizes, std::vector<std::chrono::milliseconds>(numRuns));

    toml::table config = toml::parse_file(configPath);

    const auto baseNumCellsX = config["grid"]["number_cells_x"].value_or(std::size_t{ 4 });
    const auto baseNumCellsY = config["grid"]["number_cells_y"].value_or(std::size_t{ 1 });
    const auto baseNumCellsZ = config["grid"]["number_cells_z"].value_or(std::size_t{ 1 });

    const auto scalingFactor = std::cbrt(2);

    for (auto r = 0; r < numRuns; ++r) {
        for (auto s = 0; s < numSizes; ++s) {
            const auto numCellsX = static_cast<std::size_t>(std::round(baseNumCellsX * std::pow(scalingFactor, s)));
            const auto numCellsY = static_cast<std::size_t>(std::round(baseNumCellsY * std::pow(scalingFactor, s)));
            const auto numCellsZ = static_cast<std::size_t>(std::round(baseNumCellsZ * std::pow(scalingFactor, s)));

            *(config["grid"]["number_cells_x"].as_integer()) = numCellsX;
            *(config["grid"]["number_cells_y"].as_integer()) = numCellsY;
            *(config["grid"]["number_cells_z"].as_integer()) = numCellsZ;

            ShockTube shockTube(config);

            const auto timeStart = std::chrono::high_resolution_clock::now();

            Solver solver(shockTube);
            solver.initialise();

            int numTimesteps = 0;
            while (!solver.isFinished()) {
                solver.step();
                solver.adjust();
                ++numTimesteps;
            }
            solver.finalise();

            std::cout << "Num Timesteps: " << numTimesteps << std::endl;

            const auto timeEnd = std::chrono::high_resolution_clock::now();
            const auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart);
            results[s][r] = elapsedMs;
        }
    }

    return results;
}

TEST_CASE("Shock-Tube benchmark", "[Benchmark]") {
    const auto benchmarkResults = benchmarkSolver<RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS>>(
        5, 3, "configuration/shock-tube-benchmark.toml");

    for (std::size_t s = 0; s < benchmarkResults.size(); ++s) {
        for (std::size_t r = 0; r < benchmarkResults[s].size(); ++r) {
            std::cout << "CPU: Size " << s << " run #" << r << " took " << benchmarkResults[s][r].count() << "ms"
                      << std::endl;
        }
    }
}

TEST_CASE("Shock-Tube sycl benchmark", "[Benchmark]") {
    const auto benchmarkResults = benchmarkSolver<RungeKuttaSyclSolver<ShockTube, FieldStruct, GHOST_CELLS>>(
        5, 3, "configuration/shock-tube-benchmark.toml");

    for (std::size_t s = 0; s < benchmarkResults.size(); ++s) {
        for (std::size_t r = 0; r < benchmarkResults[s].size(); ++r) {
            std::cout << "GPU: Size " << s << " run #" << r << " took " << benchmarkResults[s][r].count() << "ms"
                      << std::endl;
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

template <typename T> celerity::buffer<T, 1> uploadCelerityBuffer(const std::vector<T>& data) {
    return celerity::buffer<T, 1>(data.data(), celerity::range<1>(data.size()));
}

TEST_CASE("Celerity reduction", "[celerity]") {
    constexpr auto timedCreateQueue = []() {
        const auto _ = ScopedTimer("Creating celerity queue");
        return get_celerity_queue();
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
        constexpr auto timedBufferUpload = [](const auto& data) {
            const auto _ = ScopedTimer("Uploading buffer to GPU");
            return uploadCelerityBuffer(data);
        };

        auto buffer = timedBufferUpload(data);

        auto maxBuffer = celerity::buffer<int, 1>{ { 1 } };

        {
            const auto _ = ScopedTimer("Celerity reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
            queue.submit([=](celerity::handler& cgh) {
                auto bufferAccessor =
                    celerity::accessor{ buffer, cgh, celerity::access::one_to_one{}, celerity::read_only };

                auto maxReduction = celerity::reduction(maxBuffer, cgh, sycl::maximum<>());

                cgh.parallel_for(buffer.get_range(), maxReduction,
                                 [=](celerity::item<1> idx, auto& max) { max.combine(bufferAccessor[idx]); });
            });
        }

        queue.submit(celerity::allow_by_ref, [=, &celerityResult](celerity::handler& cgh) {
            celerity::accessor bufferAccessor{ maxBuffer, cgh, celerity::access::all{}, celerity::read_only_host_task };
            cgh.host_task(celerity::on_master_node, [=, &celerityResult] {
                celerityResult = bufferAccessor[0];
                CHECK(celerityResult == stlResult);
            });
        });

        queue.slow_full_sync();
    }
    {
        const auto _ = ScopedTimer("STL reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
        stlResult = std::accumulate(data.begin(), data.end(), std::numeric_limits<int>::lowest(), reductionOp);
    }
}
