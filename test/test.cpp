#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../src/configuration/shock-tube.h"
#include "../src/grid/grid-functions.h"
#include "../src/solver/runge-kutta-solver.h"

TEST_CASE("Shock-Tube integration test", "[IntegrationTest]") {
    ShockTube problem = ShockTube::initialiseTestProblem();

    PaddedGrid<FieldStruct, GHOST_CELLS> grid({}, problem.numberCells[Direction::DirX],
                                              problem.numberCells[Direction::DirY],
                                              problem.numberCells[Direction::DirZ]);

    RungeKuttaSolver<ShockTube> solver(grid, problem);
    solver.initialise();

    double deviationPerSteps = 1.0 + 1e-4;
    double deviationThreshold = deviationPerSteps;
    for (unsigned timeStep = 1; timeStep <= 16; timeStep++) {
        solver.step();
        if (timeStep % 4 == 0) {
            const SimpleGrid<FieldStruct> baseline =
                GridFunctions::readFromFile("test-data/step-" + std::to_string(timeStep) + ".dat");
            double averageDeviation = GridFunctions::compare(baseline, grid, false);
            REQUIRE(averageDeviation < deviationThreshold - 1.0);
            deviationThreshold *= deviationPerSteps;
        }
    }
}
