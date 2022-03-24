#include "simple-solver.h"

#include <vector>

void SimpleSolver::computeStep(Grid<double>& grid, const double deltaTime) {
    size_t xDim = grid.xDim();
    size_t yDim = grid.yDim();
    size_t zDim = grid.zDim();

    std::vector<double> temp(xDim * yDim * zDim);

    for (size_t x = 1; x < xDim - 1; x++) {
        for (size_t y = 1; y < yDim - 1; y++) {
            for (size_t z = 1; z < zDim - 1; z++) {
                temp[x * xDim * xDim + y * yDim + z] =
                    deltaTime *
                    (-6.0 * grid(x, y, z) + grid(x + 1, y, z) + grid(x - 1, y, z) +
                     grid(x, y + 1, z) + grid(x, y - 1, z) + grid(x, y, z + 1) + grid(x, y, z - 1));
            }
        }
    }

    for (size_t x = 1; x < xDim - 1; x++) {
        for (size_t y = 1; y < yDim - 1; y++) {
            for (size_t z = 1; z < zDim - 1; z++) {
                grid(x, y, z) += temp[x * xDim * xDim + y * yDim + z];
            }
        }
    }
}
