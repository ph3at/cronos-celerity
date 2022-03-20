#include <iostream>
#include <stdlib.h>

#include "grid/simple-grid.h"

using namespace std;

const size_t CELLS_PER_DIMENSION = 10;

int main(int argc, char** argv) {

    SimpleGrid grid =
        SimpleGrid::init_linear(CELLS_PER_DIMENSION, CELLS_PER_DIMENSION, CELLS_PER_DIMENSION);
    grid.print();

    return EXIT_SUCCESS;
}
