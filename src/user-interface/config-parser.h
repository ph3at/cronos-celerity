#pragma once

#include <string>
#include <toml++/toml.h>

#include "../amr/amr-parameters.h"
#include "run.h"

void parseAndRun(const std::string& configFile);

AMRParameters parseAMRConfig(const toml::table& config);

template <class ProblemType> void parseAndRun(const toml::table& config) {
    const ProblemType problem(config);
    if (config["enable_amr"].value_or(false)) {
        AMRParameters amrConfig = parseAMRConfig(config);
        runAMR(problem, amrConfig);
    } else {
        runRKS(problem);
    }
}
