#include <cstdlib>
#include <iostream>

#include "user-interface/config-parser.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Missing configuration file." << std::endl;
        return EXIT_FAILURE;
    } else {
        parseAndRun(argv[1]);
    }

    return EXIT_SUCCESS;
}
