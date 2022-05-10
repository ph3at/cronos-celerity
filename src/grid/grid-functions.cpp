#include <cmath>
#include <fstream>

#include "../misc/faces.h"
#include "grid-functions.h"

namespace GridFunctions {
bool checkNaN(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid) {
    for (unsigned x = 0; x < grid.xDim(); x++) {
        for (unsigned y = 0; y < grid.yDim(); y++) {
            for (unsigned z = 0; z < grid.zDim(); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    if (std::isnan(grid(x, y, z)[field])) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

SimpleGrid<FieldStruct> readFromFile(const std::string filename) {
    std::ifstream inputFile(filename, std::ios::in);
    if (!inputFile.is_open()) {
        std::cerr << "File \"" << filename << "\" could not be opened." << std::endl;
        return SimpleGrid<FieldStruct>({}, 0, 0, 0);
    } else {
        unsigned xDim, yDim, zDim;
        inputFile >> xDim >> yDim >> zDim;
        SimpleGrid<FieldStruct> result = SimpleGrid<FieldStruct>({}, xDim, yDim, zDim);
        double value;
        for (unsigned x = 0; x < xDim; x++) {
            for (unsigned y = 0; y < yDim; y++) {
                for (unsigned z = 0; z < zDim; z++) {
                    for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                        inputFile >> value;
                        result(x, y, z)[field] = value;
                    }
                }
            }
        }
        return result;
    }
}
}; // namespace GridFunctions
