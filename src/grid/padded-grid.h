#pragma once

#include "base-grid.h"

#include <array>
#include <vector>

template <class T> class PaddedGrid : public Grid<T> {
  public:
    PaddedGrid(const size_t xDim, const size_t yDim, const size_t zDim, size_t padding);

    T operator()(const size_t x, const size_t y, const size_t z) const;
    T& operator()(const size_t x, const size_t y, const size_t z);

  private:
    PaddedGrid();

    std::vector<T> data;
    size_t padding;
};

template <class T>
inline PaddedGrid<T>::PaddedGrid(const size_t xDim, const size_t yDim, const size_t zDim,
                                 size_t padding) {
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

template <class T>
inline T PaddedGrid<T>::operator()(const size_t x, const size_t y, const size_t z) const {
    return this->data[x * this->xSize * this->xSize + y * this->ySize + z];
}

template <class T>
inline T& PaddedGrid<T>::operator()(const size_t x, const size_t y, const size_t z) {
    return this->data[x * this->xSize * this->xSize + y * this->ySize + z];
}
