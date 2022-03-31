#include <iostream>

#include "simple-grid.h"

template <>
SimpleGrid<double> SimpleGrid<double>::initLinear(const size_t xDim, const size_t yDim,
                                                  const size_t zDim) {
    SimpleGrid grid(0.0);
    grid.xSize = xDim;
    grid.ySize = yDim;
    grid.zSize = yDim;

    grid.data.reserve(xDim * yDim * zDim);
    for (size_t x = 0; x < xDim; x++) {
        for (size_t y = 0; y < yDim; y++) {
            for (size_t z = 0; z < zDim; z++) {
                grid.data.push_back(static_cast<double>(x + y + z));
            }
        }
    }
    return grid;
}

template <> void SimpleGrid<double>::print() const {
    std::cout.width(5);
    for (size_t x = 0; x < this->xSize; x++) {
        for (size_t y = 0; y < this->ySize; y++) {
            for (size_t z = 0; z < this->zSize; z++) {
                printf(" %4.f", this->data[x * this->xSize * this->xSize + y * this->ySize + z]);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}
