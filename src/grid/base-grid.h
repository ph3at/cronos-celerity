#ifndef BASE_GRID
#define BASE_GRID

#include <cstddef>
using namespace std;

class Grid {
  public:
    virtual double* operator()(const size_t x, const size_t y, const size_t z) = 0;
    double x_dim() { return this->x_size; }
    double y_dim() { return this->y_size; }
    double z_dim() { return this->z_size; }

  protected:
    double x_size;
    double y_size;
    double z_size;
};

#endif
