#include "simple-solver.h"

#include <vector>

void SimpleSolver::solve() {
    double time = 0;
    while (time < this->timeEnd) {
        this->computeStep();
        time += this->timeDelta;
    }
}

void SimpleSolver::computeStep() {
    size_t xDim = this->grid.xDim();
    size_t yDim = this->grid.yDim();
    size_t zDim = this->grid.zDim();

    std::vector<double> temp(xDim * yDim * zDim);

    for (size_t x = 1; x < xDim - 1; x++) {
        for (size_t y = 1; y < yDim - 1; y++) {
            for (size_t z = 1; z < zDim - 1; z++) {
                temp[x * xDim * xDim + y * yDim + z] =
                    this->timeDelta *
                    (-6.0 * this->grid(x, y, z) + this->grid(x + 1, y, z) +
                     this->grid(x - 1, y, z) + this->grid(x, y + 1, z) + this->grid(x, y - 1, z) +
                     this->grid(x, y, z + 1) + this->grid(x, y, z - 1));
            }
        }
    }

    for (size_t x = 1; x < xDim - 1; x++) {
        for (size_t y = 1; y < yDim - 1; y++) {
            for (size_t z = 1; z < zDim - 1; z++) {
                this->grid(x, y, z) += temp[x * xDim * xDim + y * yDim + z];
            }
        }
    }
}
