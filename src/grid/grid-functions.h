#pragma once

#include <cassert>
#include <cmath>

#include "padded-grid.h"
#include "simple-grid.h"

#include "../configuration/constants.h"
#include "../field-wrapper/phys-fields.h"

namespace GridFunctions {
bool checkNaN(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid);
template <unsigned padding> void printGrid(const PaddedGrid<FieldStruct, padding>& grid);
template <unsigned padding1, unsigned padding2>
void compare(const PaddedGrid<FieldStruct, padding1>& baseline,
             const PaddedGrid<FieldStruct, padding2>& other);
SimpleGrid<FieldStruct> readFromFile(const std::string filename);
}; // namespace GridFunctions

template <unsigned padding>
void GridFunctions::printGrid(const PaddedGrid<FieldStruct, padding>& grid) {
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

template <unsigned padding1, unsigned padding2>
void GridFunctions::compare(const PaddedGrid<FieldStruct, padding1>& baseline,
                            const PaddedGrid<FieldStruct, padding2>& other) {
    assert(baseline.xEnd() - baseline.xStart() == other.xEnd() - other.xStart());
    assert(baseline.yEnd() - baseline.yStart() == other.yEnd() - other.yStart());
    assert(baseline.zEnd() - baseline.zStart() == other.zEnd() - other.zStart());
    const int offset = padding2 - padding1;

    double maxDiff = 0;
    double sumSquaredError = 0;
    for (int x = baseline.xStart(); x < static_cast<int>(baseline.xEnd()); x++) {
        for (int y = baseline.yStart(); y < static_cast<int>(baseline.yEnd()); y++) {
            for (int z = baseline.zStart(); z < static_cast<int>(baseline.zEnd()); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    double diff = std::abs(baseline(x, y, z)[field] -
                                           other(offset + x, offset + y, offset + z)[field]);
                    if (diff > maxDiff) {
                        maxDiff = diff;
                    }
                    sumSquaredError += diff * diff;
                }
            }
        }
    }
    double numValues = static_cast<double>(
        (baseline.xEnd() - baseline.xStart()) * (baseline.yEnd() - baseline.yStart()) *
        (baseline.zEnd() - baseline.zStart()) * NUM_PHYSICAL_FIELDS);
    double rms = std::sqrt(sumSquaredError / numValues);
    std::cout << "Root of mean squared error: " << rms << ", with a maximum deviation of "
              << maxDiff << std::endl;
}

