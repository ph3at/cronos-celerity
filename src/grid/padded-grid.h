#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "../boundary/boundary-types.h"
#include "../data-types/direction.h"
#include "../data-types/faces.h"

template <class T, unsigned padding> class PaddedGrid {
  public:
    PaddedGrid(const T defaultValue, const size_t xDim, const size_t yDim, const size_t zDim);
    PaddedGrid(const PaddedGrid<T, padding>& blueprint);

    void clear();
    void swap(PaddedGrid<T, padding>& other);

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

    T defaultValue;

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

template <class T, unsigned padding>
PaddedGrid<T, padding>::PaddedGrid(const PaddedGrid<T, padding>& blueprint)
    : defaultValue(blueprint.defaultValue), xSize(blueprint.xSize), ySize(blueprint.ySize),
      zSize(blueprint.zSize) {
    this->data.reserve(blueprint.data.size());
    for (size_t i = 0; i < blueprint.data.size(); i++) {
        this->data.push_back(blueprint.data[i]);
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
void PaddedGrid<T, padding>::swap(PaddedGrid<T, padding>& other) {
    this->data.swap(other.data);
    std::swap(this->xSize, other.xSize);
    std::swap(this->ySize, other.ySize);
    std::swap(this->zSize, other.zSize);
    this->defaultValue.swap(other.defaultValue);
}

template <class T, unsigned padding>
T& PaddedGrid<T, padding>::operator()(const size_t x, const size_t y, const size_t z) {
    return this->data[x * this->ySize * this->zSize + y * this->zSize + z];
}

template <class T, unsigned padding>
T PaddedGrid<T, padding>::operator()(const size_t x, const size_t y, const size_t z) const {
    return this->data[x * this->ySize * this->zSize + y * this->zSize + z];
}
