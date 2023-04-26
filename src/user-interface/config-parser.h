#pragma once

#include <string>
#include <toml++/toml.h>

#include "../amr/amr-parameters.h"
#include "run.h"

enum class Backend {
    CPU,
    SYCL,
    CELERITY,
};

void parseAndRun(const std::string& configFile, const Backend backend);

AMRParameters parseAMRConfig(const toml::table& config);

template <class ProblemType> void parseAndRun(const toml::table& config, const Backend backend) {
    const ProblemType problem(config);
    if (config["enable_amr"].value_or(false)) {
        AMRParameters amrConfig = parseAMRConfig(config);
        runAMR(problem, amrConfig);
    } else if (backend == Backend::CPU) {
        runRKS(problem);
    } else if (backend == Backend::SYCL) {
        runSyclRKS(problem);
    } else if (backend == Backend::CELERITY) {
        runCelerityRKS(problem);
    }
}
