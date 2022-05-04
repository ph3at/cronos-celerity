#pragma once

#include "padded-grid.h"
#include "simple-grid.h"

#include "../configuration/constants.h"
#include "../field-wrapper/phys-fields.h"

namespace GridFunctions {
bool checkNaN(PaddedGrid<FieldStruct, GHOST_CELLS>& grid);
template <unsigned padding> void printGrid(PaddedGrid<FieldStruct, padding>& grid);
}; // namespace GridFunctions

template <unsigned padding> void GridFunctions::printGrid(PaddedGrid<FieldStruct, padding>& grid) {
    std::cout.width(5);
    for (size_t x = 0; x < grid.xDim(); x++) {
        for (size_t y = 0; y < grid.yDim(); y++) {
            for (size_t z = 0; z < grid.zDim(); z++) {
                std::cout << " (";
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    printf("%8.4f", grid(x, y, z)[field]);
                    if (field < NUM_PHYSICAL_FIELDS - 1) {
                        std::cout << ",";
                    }
                }
                std::cout << ")";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}
