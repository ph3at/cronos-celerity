#include "simple-solver.h"

#include <vector>

void SimpleSolver::compute_step(Grid* grid, const double time_step) {
    size_t x_dim = grid->x_dim();
    size_t y_dim = grid->y_dim();
    size_t z_dim = grid->z_dim();

    vector<double> temp(x_dim * y_dim * z_dim);

    for (size_t x = 1; x < x_dim - 1; x++) {
        for (size_t y = 1; y < y_dim - 1; y++) {
            for (size_t z = 1; z < z_dim - 1; z++) {
                temp[x * x_dim * x_dim + y * y_dim + z] =
                    time_step *
                    (-6.0 * *(*grid)(x, y, z) + *(*grid)(x + 1, y, z) + *(*grid)(x - 1, y, z) +
                     *(*grid)(x, y + 1, z) + *(*grid)(x, y - 1, z) + *(*grid)(x, y, z + 1) +
                     *(*grid)(x, y, z - 1));
            }
        }
    }

    for (size_t x = 1; x < x_dim - 1; x++) {
        for (size_t y = 1; y < y_dim - 1; y++) {
            for (size_t z = 1; z < z_dim - 1; z++) {
                *(*grid)(x, y, z) += temp[x * x_dim * x_dim + y * y_dim + z];
            }
        }
    }
}
