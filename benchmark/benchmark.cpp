#include <toml++/toml.h>

#include <celerity.h>
#include <sycl/sycl.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "../src/configuration/shock-tube.h"
#include "../src/solver/runge-kutta-celerity-solver.h"
#include "../src/solver/runge-kutta-solver.h"
#include "../src/solver/runge-kutta-sycl-solver.h"

template <typename Solver, typename... Args>
auto benchmarkSolver(const int numSizes, const int numRuns, const std::string_view configPath, Args&... args) {
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

            Solver solver(shockTube, args...);
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

void cpuBenchmark() {
    const auto benchmarkResults = benchmarkSolver<RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS>>(
        20, 10, "configuration/shock-tube-benchmark.toml");

    for (std::size_t s = 0; s < benchmarkResults.size(); ++s) {
        for (std::size_t r = 0; r < benchmarkResults[s].size(); ++r) {
            std::cout << "CPU: Size " << s << " run #" << r << " took " << benchmarkResults[s][r].count() << "ms"
                      << std::endl;
        }
    }
}

void syclBenchmark() {
    const auto benchmarkResults = benchmarkSolver<RungeKuttaSyclSolver<ShockTube, FieldStruct, GHOST_CELLS>>(
        24, 10, "configuration/shock-tube-benchmark.toml");

    for (std::size_t s = 0; s < benchmarkResults.size(); ++s) {
        for (std::size_t r = 0; r < benchmarkResults[s].size(); ++r) {
            std::cout << "Sycl GPU: Size " << s << " run #" << r << " took " << benchmarkResults[s][r].count() << "ms"
                      << std::endl;
        }
    }
}

void celerityBenchmark() {
    auto queue = celerity::distr_queue();
    const auto benchmarkResults = benchmarkSolver<RungeKuttaCeleritySolver<ShockTube, FieldStruct, GHOST_CELLS>>(
        24, 10, "configuration/shock-tube-benchmark.toml", queue);

    queue.submit(celerity::allow_by_ref, [=, &benchmarkResults](celerity::handler& cgh) {
        cgh.host_task(celerity::on_master_node, [=, &benchmarkResults] {
            for (std::size_t s = 0; s < benchmarkResults.size(); ++s) {
                for (std::size_t r = 0; r < benchmarkResults[s].size(); ++r) {
                    CELERITY_INFO("Celerity GPU: Size {} run #{} took {}ms", s, r, benchmarkResults[s][r].count());
                }
            }
        });
    });
    queue.slow_full_sync();
}

void printUsage(const std::string_view name) {
    std::cout << "Usage: " << name << " cpu|sycl|celerity" << std::endl;
}

int main(int argc, char* argv[]) {
    const auto arguments = std::vector<std::string>{ argv, argv + argc };

    if (arguments.size() != 2) {
        printUsage(arguments[0]);
        return EXIT_FAILURE;
    }

    const auto benchmark = arguments[1];

    if (benchmark == "cpu") {
        cpuBenchmark();
    } else if (benchmark == "sycl") {
        syclBenchmark();
    } else if (benchmark == "celerity") {
        celerityBenchmark();
    } else {
        std::cout << "ERROR - Unknown benchmark type '" << benchmark << "'" << std::endl;
        printUsage(arguments[1]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
