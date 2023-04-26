#include <iostream>

#include "../configuration/shock-tube.h"
#include "config-parser.h"
#include "value-parser.h"

void parseAndRun(const std::string& configFile, const Backend backend) {
    try {
        toml::table config = toml::parse_file(configFile);
        const std::optional<std::string> maybeProblemType = config["problem"].value<std::string>();
        if (maybeProblemType.has_value()) {
            const std::string& problemType = maybeProblemType.value();
            if (problemType.compare("shock tube") == 0) {
                std::cout << "--------------- Running shock tube ---------------" << std::endl;
                parseAndRun<ShockTube>(config, backend);
            } else {
                std::cerr << "Unknown problem type: " << problemType << std::endl;
                exit(4);
            }
        } else {
            std::cerr << "Missing problem type." << std::endl;
            exit(3);
        }
    } catch (const toml::parse_error& err) {
        std::cerr << "Unable to parse file " << configFile << ":" << std::endl << err << std::endl;
        exit(2);
    }
}

AMRParameters parseAMRConfig(const toml::table& config) {
    std::cout << std::endl << "AMR parameters:" << std::endl;
    const toml::node_view<const toml::node>& amr = config["amr"];
    const AMRParameters amrConfig = { parseValue<unsigned>(amr, "refinement_factor", 4),
                                      parseValue<unsigned>(amr, "refinement_interval", 4),
                                      parseValue<unsigned>(amr, "buffer_size", 4),
                                      parseValue<double>(amr, "efficiency_threshold", 0.6),
                                      parseValue<double>(amr, "truncation_error_threshold", 1e-4),
                                      parseValue<unsigned>(amr, "max_refinement_depth", 1) };
    return amrConfig;
}
