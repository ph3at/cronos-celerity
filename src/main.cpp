#include "user-interface/run.h"
#include <cstdlib>

int main(int argc, char** argv) {
    if (argc < 2 || std::string("-amrOff").compare(argv[1]) != 0) {
        runAMR<ShockTube>();
    } else {
        runRK<ShockTube>();
    }

    return EXIT_SUCCESS;
}
