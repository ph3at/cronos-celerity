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
auto benchmarkSolver(const int size, const std::string_view configPath, Args&... args) {
    toml::table config = toml::parse_file(configPath);

    const auto baseNumCellsX = config["grid"]["number_cells_x"].value_or(std::size_t{ 4 });
    const auto baseNumCellsY = config["grid"]["number_cells_y"].value_or(std::size_t{ 1 });
    const auto baseNumCellsZ = config["grid"]["number_cells_z"].value_or(std::size_t{ 1 });

    const auto scalingFactor = std::cbrt(2);

    const auto numCellsX = static_cast<std::size_t>(std::round(baseNumCellsX * std::pow(scalingFactor, size)));
    const auto numCellsY = static_cast<std::size_t>(std::round(baseNumCellsY * std::pow(scalingFactor, size)));
    const auto numCellsZ = static_cast<std::size_t>(std::round(baseNumCellsZ * std::pow(scalingFactor, size)));

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

    const auto timeEnd = std::chrono::high_resolution_clock::now();
    const auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart);

    std::cout << "Num Timesteps: " << numTimesteps << std::endl;

    return elapsedMs;
}

void cpuBenchmark(const int size) {
    const auto benchmarkTimeMs = benchmarkSolver<RungeKuttaSolver<ShockTube, FieldStruct, GHOST_CELLS>>(
        size, "configuration/shock-tube-benchmark.toml");

    std::cerr << "arch,size,time" << std::endl;
    std::cerr << "cpu," << size << "," << benchmarkTimeMs.count() << std::endl;
}

void syclBenchmark(const int size) {
    const auto benchmarkTimeMs = benchmarkSolver<RungeKuttaSyclSolver<ShockTube, FieldStruct, GHOST_CELLS>>(
        size, "configuration/shock-tube-benchmark.toml");

    std::cerr << "arch,size,time" << std::endl;
    std::cerr << "sycl," << size << "," << benchmarkTimeMs.count() << std::endl;
}

void celerityBenchmark(const int size) {
    auto queue = celerity::distr_queue();
    const auto benchmarkTimeMs = benchmarkSolver<RungeKuttaCeleritySolver<ShockTube, FieldStruct, GHOST_CELLS>>(
        size, "configuration/shock-tube-benchmark.toml", queue);

    queue.submit(celerity::allow_by_ref, [=, &benchmarkTimeMs](celerity::handler& cgh) {
        cgh.host_task(celerity::on_master_node, [=, &benchmarkTimeMs] {
            std::cerr << "arch,size,time" << std::endl;
            std::cerr << "celerity," << size << "," << benchmarkTimeMs.count() << std::endl;
        });
    });
    queue.slow_full_sync();
}

void printUsage(const std::string_view name) {
    std::cout << "Usage: " << name << " cpu|sycl|celerity size" << std::endl;
}

int main(int argc, char* argv[]) {
    const auto arguments = std::vector<std::string>{ argv, argv + argc };

    if (arguments.size() != 3) {
        printUsage(arguments[0]);
        return EXIT_FAILURE;
    }

    const auto benchmark = arguments[1];
    const auto size = std::stoi(arguments[2]);

    if (size < 0 || size > 100) {
        std::cout << "ERROR - Invalid size value '" << size << "' should be between [0-100]" << std::endl;
        printUsage(arguments[0]);
        return EXIT_FAILURE;
    }

    if (benchmark == "cpu") {
        cpuBenchmark(size);
    } else if (benchmark == "sycl") {
        syclBenchmark(size);
    } else if (benchmark == "celerity") {
        celerityBenchmark(size);
    } else {
        std::cout << "ERROR - Unknown benchmark type '" << benchmark << "'" << std::endl;
        printUsage(arguments[0]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
