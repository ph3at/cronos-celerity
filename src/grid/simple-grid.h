#ifndef GRID
#define GRID

#include <vector>

using namespace std;

typedef vector<double> Grid;

Grid init_linear_grid(size_t x_dim, size_t y_dim, size_t z_dim);

void print_grid(Grid& grid, size_t x_dim, size_t y_dim, size_t z_dim);

#endif
