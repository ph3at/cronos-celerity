#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <toml++/toml.h>

#include <CL/sycl.hpp>

#include <algorithm>
#include <bit>
#include <chrono>
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
    for (unsigned timeStep = 1; timeStep <= 16; timeStep++) {
        solver.step();
        syclSolver.step();
        solver.adjust();
        syclSolver.adjust();
        const auto& baseline = solver.grid;
        double averageDeviation = GridFunctions::compare(baseline, syclSolver.grid(), false, false);
        CAPTURE(timeStep);
        CHECK(averageDeviation == 0);
        REQUIRE(averageDeviation < deviationThreshold - 1.0);
        deviationThreshold *= deviationPerStep;
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
        auto queue = cl::sycl::queue();
        auto itemRange = cl::sycl::range<3>(NUM_ELEMENTS[0], NUM_ELEMENTS[1], NUM_ELEMENTS[2]);
        auto totalRange = cl::sycl::range<3>(DATA_SIZE[0], DATA_SIZE[1], DATA_SIZE[2]);
        auto bufferA = cl::sycl::buffer<int, 3>(dataA.data(), totalRange);
        auto bufferB = cl::sycl::buffer<int, 3>(dataB.data(), totalRange);
        auto bufferC = cl::sycl::buffer<int, 3>(dataC.data(), totalRange);

        queue.submit([&](cl::sycl::handler& cgh) {
            auto accA = bufferA.template get_access<cl::sycl::access::mode::read>(cgh);
            auto accB = bufferB.template get_access<cl::sycl::access::mode::read>(cgh);
            auto accC = bufferC.template get_access<cl::sycl::access::mode::discard_write>(cgh);

            cgh.parallel_for(itemRange, [=](const cl::sycl::id<3> id) {
                const auto offset = cl::sycl::id<3>(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);
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

template <typename T> cl::sycl::buffer<T, 1> uploadBuffer(const std::vector<T>& data) {
    return cl::sycl::buffer<T, 1>(data.data(), cl::sycl::range<1>(data.size()));
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

TEST_CASE("Sycl reduction", "[sycl]") {
    constexpr auto timedCreateQueue = []() {
        const auto _ = ScopedTimer("Creating SYCL queue");
        return cl::sycl::queue([](const cl::sycl::exception_list exceptions) {
            try {
                for (const auto& e : exceptions) {
                    std::rethrow_exception(e);
                }
            } catch (const cl::sycl::exception& e) {
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
            syclResult = utils::sycl::reduce(queue, buffer, reductionOp, std::numeric_limits<int>::lowest());
        }
    }
    {
        const auto _ = ScopedTimer("STL reduction of " + std::to_string(NUM_ELEMENTS) + " elements");
        stlResult = std::accumulate(data.begin(), data.end(), std::numeric_limits<int>::lowest(), reductionOp);
    }

    CHECK(syclResult == stlResult);
    std::cout << std::endl;
}
