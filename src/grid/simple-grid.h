#pragma once

#include <iostream>
#include <vector>

#include "base-grid.h"

template <class T> class SimpleGrid : public Grid<T> {
  public:
    SimpleGrid(const T defaultValue, const size_t xDim, const size_t yDim, const size_t zDim);
    static SimpleGrid<double> initLinear(const size_t xDim, const size_t yDim, const size_t zDim);

    T operator()(const size_t x, const size_t y, const size_t z) const;
    T& operator()(const size_t x, const size_t y, const size_t z);

    void clear();

    void print() const;
    void setBorderConst(const T borderValue);

    size_t xStart() const { return 0; }
    size_t yStart() const { return 0; }
    size_t zStart() const { return 0; }
    size_t xEnd() const { return this->xSize; }
    size_t yEnd() const { return this->ySize; }
    size_t zEnd() const { return this->zSize; }

  private:
    SimpleGrid(const T defaultValue);
    std::vector<T> data;
};

template <class T> inline SimpleGrid<T>::SimpleGrid(const T defaultValue) : Grid<T>(defaultValue) {}

template <class T>
inline SimpleGrid<T>::SimpleGrid(const T defaultValue, const size_t xDim, const size_t yDim,
                                 const size_t zDim)
    : Grid<T>(defaultValue) {
    this->xSize = xDim;
    this->ySize = yDim;
    this->zSize = zDim;

    const size_t size = this->xSize * this->ySize * this->zSize;
    this->data.reserve(size);
    for (size_t _ = 0; _ < size; _++) {
        this->data.push_back(this->defaultValue);
    }
}

template <class T> void SimpleGrid<T>::clear() {
    for (unsigned x = this->xStart(); x < this->xEnd(); x++) {
        for (unsigned y = this->yStart(); y < this->yEnd(); y++) {
            for (unsigned z = this->zStart(); z < this->zEnd(); z++) {
                this->data[x * this->xSize * this->xSize + y * this->ySize + z] =
                    this->defaultValue;
            }
        }
    }
}

template <class T>
inline T SimpleGrid<T>::operator()(const size_t x, const size_t y, const size_t z) const {
    return this->data[x * this->xSize * this->xSize + y * this->ySize + z];
}

template <class T>
inline T& SimpleGrid<T>::operator()(const size_t x, const size_t y, const size_t z) {
    return this->data[x * this->xSize * this->xSize + y * this->ySize + z];
}

template <class T> inline void SimpleGrid<T>::print() const {
    std::cout.width(5);
    for (size_t x = 0; x < this->xSize; x++) {
        for (size_t y = 0; y < this->ySize; y++) {
            for (size_t z = 0; z < this->zSize; z++) {
                std::cout << " " << this->data[x * this->xSize * this->xSize + y * this->ySize + z];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}

template <class T> inline void SimpleGrid<T>::setBorderConst(const T borderValue) {
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
