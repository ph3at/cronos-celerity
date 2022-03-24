#pragma once

#include <cstddef>

class Grid {
  public:
    virtual double operator()(const size_t x, const size_t y, const size_t z) const = 0;
    virtual double& operator()(const size_t x, const size_t y, const size_t z) = 0;
    double xDim() const { return this->xSize; }
    double yDim() const { return this->ySize; }
    double zDim() const { return this->zSize; }

  protected:
    double xSize;
    double ySize;
    double zSize;
};
