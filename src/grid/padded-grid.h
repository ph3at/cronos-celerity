#pragma once

#include "base-grid.h"

#include <array>
#include <vector>

template <class T> class PaddedGrid : public Grid<T> {
  public:
    PaddedGrid(const T defaultValue, const size_t xDim, const size_t yDim, const size_t zDim,
               size_t padding);

    void clear();

    T operator()(const size_t x, const size_t y, const size_t z) const;
    T& operator()(const size_t x, const size_t y, const size_t z);

    size_t xStart() const { return padding; }
    size_t yStart() const { return padding; }
    size_t zStart() const { return padding; }
    size_t xEnd() const { return this->xSize - padding; }
    size_t yEnd() const { return this->ySize - padding; }
    size_t zEnd() const { return this->zSize - padding; }

  private:
    PaddedGrid();

    std::vector<T> data;
    size_t padding;
};

template <class T>
PaddedGrid<T>::PaddedGrid(const T defaultValue, const size_t xDim, const size_t yDim,
                          const size_t zDim, size_t padding)
    : Grid<T>(defaultValue) {
    this->xSize = xDim + 2 * padding;
    this->ySize = yDim + 2 * padding;
    this->zSize = zDim + 2 * padding;

    this->padding = padding;

    const size_t size = this->xSize * this->ySize * this->zSize;
    this->data.reserve(size);
    for (size_t _ = 0; _ < size; _++) {
        this->data.push_back(T());
    }
}

template <class T> void PaddedGrid<T>::clear() {
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
T PaddedGrid<T>::operator()(const size_t x, const size_t y, const size_t z) const {
    return this->data[x * this->xSize * this->xSize + y * this->ySize + z];
}

template <class T> T& PaddedGrid<T>::operator()(const size_t x, const size_t y, const size_t z) {
    return this->data[x * this->xSize * this->xSize + y * this->ySize + z];
}
