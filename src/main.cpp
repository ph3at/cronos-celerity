#include <iostream>
#include <stdlib.h>

#include "grid/simple-grid.h"

using namespace std;

const size_t CELLS_PER_DIMENSION = 10;

int main(int argc, char** argv) {

    Grid grid = init_linear_grid(CELLS_PER_DIMENSION, CELLS_PER_DIMENSION, CELLS_PER_DIMENSION);

    print_grid(grid, CELLS_PER_DIMENSION, CELLS_PER_DIMENSION, CELLS_PER_DIMENSION);

    return EXIT_SUCCESS;
}
