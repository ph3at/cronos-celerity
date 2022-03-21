#include <iostream>

#include "simple-grid.h"

SimpleGrid::SimpleGrid(){};

SimpleGrid::SimpleGrid(const size_t x_dim, const size_t y_dim, const size_t z_dim) {
    this->x_size = x_dim;
    this->y_size = y_dim;
    this->z_size = z_dim;

    this->data.assign(x_dim * y_dim * z_dim, 0.0);
}

SimpleGrid SimpleGrid::init_linear(const size_t x_dim, const size_t y_dim, const size_t z_dim) {
    SimpleGrid grid;
    grid.x_size = x_dim;
    grid.y_size = y_dim;
    grid.z_size = y_dim;

    grid.data.reserve(x_dim * y_dim * z_dim);
    for (size_t x = 0; x < x_dim; x++) {
        for (size_t y = 0; y < y_dim; y++) {
            for (size_t z = 0; z < z_dim; z++) {
                grid.data.push_back(double(x + y + z));
            }
        }
    }
    return grid;
}

void SimpleGrid::print() {
    cout.width(5);
    for (size_t x = 0; x < this->x_size; x++) {
        for (size_t y = 0; y < this->y_size; y++) {
            for (size_t z = 0; z < this->z_size; z++) {
                printf(" %4.f", this->data[x * this->x_size * this->x_size + y * this->y_size + z]);
            }
            cout << endl;
        }
        cout << endl << endl;
    }
}

double& SimpleGrid::operator()(const size_t x, const size_t y, const size_t z) {
    return this->data[x * this->x_size * this->x_size + y * this->y_size + z];
}

void SimpleGrid::set_border_const(const double border_value) {
    for (size_t x = 0; x < this->x_size; x++) {
        for (size_t y = 0; y < this->y_size; y++) {
            this->data[x * this->x_size * this->x_size + y * this->y_size] = border_value;
            this->data[x * this->x_size * this->x_size + y * this->y_size + this->z_size - 1] =
                border_value;
        }
    }

    for (size_t x = 0; x < this->x_size; x++) {
        for (size_t z = 0; z < this->z_size; z++) {
            this->data[x * this->x_size * this->x_size + z] = border_value;
            this->data[x * this->x_size * this->x_size + (this->y_size - 1) * this->y_size + z] =
                border_value;
        }
    }

    for (size_t y = 0; y < this->y_size; y++) {
        for (size_t z = 0; z < this->z_size; z++) {
            this->data[y * this->y_size + z] = border_value;
            this->data[(this->x_size - 1) * this->x_size * this->x_size + y * this->y_size + z] =
                border_value;
        }
    }
}
