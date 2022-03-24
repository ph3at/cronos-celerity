#include <iostream>

#include "simple-grid.h"

SimpleGrid::SimpleGrid(){};

SimpleGrid::SimpleGrid(const size_t xDim, const size_t yDim, const size_t zDim) {
    this->xSize = xDim;
    this->ySize = yDim;
    this->zSize = zDim;

    this->data.assign(xDim * yDim * zDim, 0.0);
}

SimpleGrid SimpleGrid::initLinear(const size_t xDim, const size_t yDim, const size_t zDim) {
    SimpleGrid grid;
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

double SimpleGrid::operator()(const size_t x, const size_t y, const size_t z) const {
    return this->data[x * this->xSize * this->xSize + y * this->ySize + z];
}

double& SimpleGrid::operator()(const size_t x, const size_t y, const size_t z) {
    return this->data[x * this->xSize * this->xSize + y * this->ySize + z];
}

void SimpleGrid::print() const {
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

void SimpleGrid::setBorderConst(const double borderValue) {
    for (size_t x = 0; x < this->xSize; x++) {
        for (size_t y = 0; y < this->ySize; y++) {
            this->data[x * this->xSize * this->xSize + y * this->ySize] = borderValue;
            this->data[x * this->xSize * this->xSize + y * this->ySize + this->zSize - 1] =
                borderValue;
        }
    }

    for (size_t x = 0; x < this->xSize; x++) {
        for (size_t z = 0; z < this->zSize; z++) {
            this->data[x * this->xSize * this->xSize + z] = borderValue;
            this->data[x * this->xSize * this->xSize + (this->ySize - 1) * this->ySize + z] =
                borderValue;
        }
    }

    for (size_t y = 0; y < this->ySize; y++) {
        for (size_t z = 0; z < this->zSize; z++) {
            this->data[y * this->ySize + z] = borderValue;
            this->data[(this->xSize - 1) * this->xSize * this->xSize + y * this->ySize + z] =
                borderValue;
        }
    }
}
