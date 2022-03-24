#pragma once

#include <cstddef>

template <class T> class Grid {
  public:
    virtual T operator()(const size_t x, const size_t y, const size_t z) const = 0;
    virtual T& operator()(const size_t x, const size_t y, const size_t z) = 0;
    size_t xDim() const { return this->xSize; }
    size_t yDim() const { return this->ySize; }
    size_t zDim() const { return this->zSize; }

  protected:
    size_t xSize;
    size_t ySize;
    size_t zSize;
};
