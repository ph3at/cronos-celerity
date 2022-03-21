#ifndef GRID
#define GRID

#include <vector>

#include "base-grid.h"

using namespace std;

class SimpleGrid : public Grid {
  public:
    SimpleGrid(const size_t x_dim, const size_t y_dim, const size_t z_dim);
    static SimpleGrid init_linear(const size_t x_dim, const size_t y_dim, const size_t z_dim);
    void print();
    double& operator()(const size_t x, const size_t y, const size_t z);
    void set_border_const(const double border_value);

  private:
    SimpleGrid();
    vector<double> data;
};

#endif
