#pragma once

#include <cassert>
#include <cmath>

#include "padded-grid.h"
#include "simple-grid.h"

#include "../configuration/constants.h"
#include "../data-types/phys-fields.h"

namespace GridFunctions {
bool checkNaN(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid);
template <unsigned padding> void printGrid(const PaddedGrid<FieldStruct, padding>& grid);
template <unsigned padding1, unsigned padding2>
double compare(const PaddedGrid<FieldStruct, padding1>& baseline,
               const PaddedGrid<FieldStruct, padding2>& other, const bool doOutput);
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
double GridFunctions::compare(const PaddedGrid<FieldStruct, padding1>& baseline,
                              const PaddedGrid<FieldStruct, padding2>& other, const bool doOutput) {
    assert(baseline.xEnd() - baseline.xStart() == other.xEnd() - other.xStart());
    assert(baseline.yEnd() - baseline.yStart() == other.yEnd() - other.yStart());
    assert(baseline.zEnd() - baseline.zStart() == other.zEnd() - other.zStart());
    const int offset = padding2 - padding1;

    double maxDiff = 0.0;
    double sumSquaredError = 0.0;
    double averageDeviation = 0.0;
    for (int x = baseline.xStart(); x < static_cast<int>(baseline.xEnd()); x++) {
        for (int y = baseline.yStart(); y < static_cast<int>(baseline.yEnd()); y++) {
            for (int z = baseline.zStart(); z < static_cast<int>(baseline.zEnd()); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    double correct = baseline(x, y, z)[field];
                    double actual = other(offset + x, offset + y, offset + z)[field];
                    double diff = std::abs(correct - actual);
                    if (diff > maxDiff) {
                        maxDiff = diff;
                    }
                    if (diff != 0) {
                        double factor;
                        if (correct != 0) {
                            factor = correct;
                        } else {
                            factor = actual;
                        }
                        double deviation = diff / factor;
                        if (doOutput && deviation > 0.05) {
                            std::cerr << "(" << x << "," << y << "," << z << ")[" << field << "]: ("
                                      << correct << "," << actual << ")" << std::endl;
                        }
                        averageDeviation += deviation;
                    }
                    sumSquaredError += diff * diff;
                }
            }
        }
    }
    double numValues = static_cast<double>(
        (baseline.xEnd() - baseline.xStart()) * (baseline.yEnd() - baseline.yStart()) *
        (baseline.zEnd() - baseline.zStart()) * NUM_PHYSICAL_FIELDS);
    averageDeviation /= numValues;
    if (doOutput) {
        double rms = std::sqrt(sumSquaredError / numValues);
        std::cout << "Root of mean squared error: " << rms << ", with a maximum deviation of "
                  << maxDiff << std::endl
                  << " and average percentile deviation " << averageDeviation * 100.0 << ".";
    }
    return averageDeviation;
}

