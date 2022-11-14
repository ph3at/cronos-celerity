#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <toml++/toml.h>

#include <CL/sycl.hpp>

#include "../src/configuration/shock-tube.h"
#include "../src/grid/grid-functions.h"
#include "../src/solver/runge-kutta-solver.h"
#include "../src/solver/runge-kutta-sycl-solver.h"

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
    constexpr auto DATA_SIZE = 1000;

    auto dataA = std::vector<float>(DATA_SIZE);
    auto dataB = std::vector<float>(DATA_SIZE);
    auto dataC = std::vector<float>(DATA_SIZE, -1);

    for (std::size_t i = 0; i < DATA_SIZE; ++i) {
        const auto value = static_cast<float>(i) / 100;
        dataA[value];
        dataB[-value];
    }

    for (const auto& elem : dataC) {
        CHECK(elem == -1);
    }

    {
        auto queue = cl::sycl::queue();
        auto itemRange = cl::sycl::range<1>(DATA_SIZE);
        auto bufferA = cl::sycl::buffer<float, 1>(dataA.data(), itemRange);
        auto bufferB = cl::sycl::buffer<float, 1>(dataB.data(), itemRange);
        auto bufferC = cl::sycl::buffer<float, 1>(dataC.data(), itemRange);

        queue.submit([&](cl::sycl::handler& cgh) {
            auto accA = bufferA.template get_access<cl::sycl::access::mode::read>(cgh);
            auto accB = bufferB.template get_access<cl::sycl::access::mode::read>(cgh);
            auto accC = bufferC.template get_access<cl::sycl::access::mode::discard_write>(cgh);

            cgh.parallel_for(itemRange,
                             [=](cl::sycl::id<1> id) { accC[id] = accA[id] + accB[id]; });
        });
    }

    for (const auto& elem : dataC) {
        CHECK(elem == 0);
    }
}
