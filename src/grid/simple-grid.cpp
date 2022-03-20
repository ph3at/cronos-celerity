#include <iostream>

#include "simple-grid.h"

Grid init_linear_grid(size_t x_dim, size_t y_dim, size_t z_dim) {
    Grid grid;
    grid.reserve(x_dim * y_dim * z_dim);
    for (size_t x = 0; x < x_dim; x++) {
        for (size_t y = 0; y < y_dim; y++) {
            for (size_t z = 0; z < z_dim; z++) {
                grid.push_back(double(x + y + z));
            }
        }
    }
    return grid;
}

void print_grid(Grid& grid, size_t x_dim, size_t y_dim, size_t z_dim) {
    cout.width(5);
    for (size_t x = 0; x < x_dim; x++) {
        for (size_t y = 0; y < y_dim; y++) {
            for (size_t z = 0; z < z_dim; z++) {
                printf(" %4.f", grid[x * x_dim * x_dim + y * y_dim + z]);
            }
            cout << endl;
        }
        cout << endl << endl;
    }
}
