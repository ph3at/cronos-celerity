#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../src/configuration/shock-tube.h"
#include "../src/grid/grid-functions.h"
#include "../src/solver/runge-kutta-solver.h"

TEST_CASE("Shock-Tube integration test", "[IntegrationTest]") {
    std::pair<PaddedGrid<FieldStruct, GHOST_CELLS>, ShockTube> shockTube =
        ShockTube::initialiseTestProblem();

    RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS> solver(shockTube.first, shockTube.second);
    solver.initialise();

    double deviationPerStep = 1.0005;
    double deviationThreshold = deviationPerStep;
    for (unsigned timeStep = 1; timeStep <= 16; timeStep++) {
        solver.step();
        solver.adjust();
        const SimpleGrid<FieldStruct> baseline =
            GridFunctions::readFromFile("test-data/step-" + std::to_string(timeStep) + ".dat");
        double averageDeviation = GridFunctions::compare(baseline, shockTube.first, false, false);
        REQUIRE(averageDeviation < deviationThreshold - 1.0);
        deviationThreshold *= deviationPerStep;
    }
}
