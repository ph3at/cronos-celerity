#ifndef GRID
#define GRID

#include <vector>

#include "base-grid.h"

class SimpleGrid : public Grid {
  public:
    SimpleGrid(const size_t xDim, const size_t yDim, const size_t zDim);
    static SimpleGrid initLinear(const size_t xDim, const size_t yDim, const size_t zDim);

    double operator()(const size_t x, const size_t y, const size_t z) const;
    double& operator()(const size_t x, const size_t y, const size_t z);

    void print() const;
    void setBorderConst(const double borderValue);

  private:
    SimpleGrid();
    std::vector<double> data;
};

#endif
