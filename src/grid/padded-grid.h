#pragma once

#include <array>
#include <cstddef>
#include <vector>

template <class T, unsigned padding> class PaddedGrid {
  public:
    PaddedGrid(const T defaultValue, const size_t xDim, const size_t yDim, const size_t zDim);

    void clear();

    T& operator()(const size_t x, const size_t y, const size_t z);
    T operator()(const size_t x, const size_t y, const size_t z) const;

    size_t xDim() const { return this->xSize; }
    size_t yDim() const { return this->ySize; }
    size_t zDim() const { return this->zSize; }
    size_t xStart() const { return padding; }
    size_t yStart() const { return padding; }
    size_t zStart() const { return padding; }
    size_t xEnd() const { return this->xSize - padding; }
    size_t yEnd() const { return this->ySize - padding; }
    size_t zEnd() const { return this->zSize - padding; }

    const T defaultValue;

  private:
    PaddedGrid();

    std::vector<T> data;
    size_t xSize;
    size_t ySize;
    size_t zSize;
};

template <class T, unsigned padding>
PaddedGrid<T, padding>::PaddedGrid(const T defaultValue, const size_t xDim, const size_t yDim,
                                   const size_t zDim)
    : defaultValue(defaultValue) {
    this->xSize = xDim + 2 * padding;
    this->ySize = yDim + 2 * padding;
    this->zSize = zDim + 2 * padding;

    const size_t size = this->xSize * this->ySize * this->zSize;
    this->data.reserve(size);
    for (size_t _ = 0; _ < size; _++) {
        this->data.push_back(T(defaultValue));
    }
}

template <class T, unsigned padding> void PaddedGrid<T, padding>::clear() {
    for (unsigned x = this->xStart(); x < this->xEnd(); x++) {
        for (unsigned y = this->yStart(); y < this->yEnd(); y++) {
            for (unsigned z = this->zStart(); z < this->zEnd(); z++) {
                this->data[x * this->ySize * this->zSize + y * this->zSize + z] =
                    this->defaultValue;
            }
        }
    }
}

template <class T, unsigned padding>
T& PaddedGrid<T, padding>::operator()(const size_t x, const size_t y, const size_t z) {
    return this->data[x * this->ySize * this->zSize + y * this->zSize + z];
}

template <class T, unsigned padding>
T PaddedGrid<T, padding>::operator()(const size_t x, const size_t y, const size_t z) const {
    return this->data[x * this->ySize * this->zSize + y * this->zSize + z];
}
