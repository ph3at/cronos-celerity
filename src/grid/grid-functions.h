#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>

#include <celerity.h>
#include <sycl/sycl.hpp>

#include "padded-grid.h"
#include "simple-grid.h"

#include "../configuration/constants.h"
#include "../data-types/phys-fields.h"

namespace GridFunctions {
bool checkNaN(const PaddedGrid<FieldStruct, GHOST_CELLS>& grid);
template <unsigned padding> void printGrid(const PaddedGrid<FieldStruct, padding>& grid);
template <unsigned padding1, unsigned padding2>
double compare(const PaddedGrid<FieldStruct, padding1>& baseline, const PaddedGrid<FieldStruct, padding2>& other,
               const bool doOutput, const bool verbose);
SimpleGrid<FieldStruct> readFromFile(const std::string filename);
template <unsigned padding> void writeToFile(const std::string fileName, const PaddedGrid<FieldStruct, padding>& grid);
}; // namespace GridFunctions

template <unsigned padding> void GridFunctions::printGrid(const PaddedGrid<FieldStruct, padding>& grid) {
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
                              const PaddedGrid<FieldStruct, padding2>& other, const bool doOutput, const bool verbose) {
    assert(baseline.xEnd() - baseline.xStart() == other.xEnd() - other.xStart());
    assert(baseline.yEnd() - baseline.yStart() == other.yEnd() - other.yStart());
    assert(baseline.zEnd() - baseline.zStart() == other.zEnd() - other.zStart());
    const int offset = padding2 - padding1;

    double maxDiff = 0.0;
    double sumSquaredError = 0.0;
    double averageDeviation = 0.0;
    int maxX = 0, maxY = 0, maxZ = 0, fieldMax = 0;
    for (int x = baseline.xStart(); x < static_cast<int>(baseline.xEnd()); x++) {
        for (int y = baseline.yStart(); y < static_cast<int>(baseline.yEnd()); y++) {
            for (int z = baseline.zStart(); z < static_cast<int>(baseline.zEnd()); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    double correct = baseline(x, y, z)[field];
                    double actual = other(offset + x, offset + y, offset + z)[field];
                    double diff = std::abs(correct - actual);
                    if (diff > maxDiff) {
                        maxDiff = diff;
                        maxX = x;
                        maxY = y;
                        maxZ = z;
                        fieldMax = field;
                    }
                    if (diff != 0) {
                        double factor;
                        if (correct != 0) {
                            factor = correct;
                        } else {
                            factor = actual;
                        }
                        double deviation = diff / std::abs(factor);
                        if (doOutput && verbose && deviation > 0.01) {
                            std::cerr << "(" << x << "," << y << "," << z << ")[" << field << "]: (" << correct << ","
                                      << actual << ")" << std::endl;
                        }
                        averageDeviation += deviation;
                    }
                    sumSquaredError += diff * diff;
                }
            }
        }
    }
    double numValues =
        static_cast<double>((baseline.xEnd() - baseline.xStart()) * (baseline.yEnd() - baseline.yStart()) *
                            (baseline.zEnd() - baseline.zStart()) * NUM_PHYSICAL_FIELDS);
    averageDeviation /= numValues;
    if (doOutput) {
        double rms = std::sqrt(sumSquaredError / numValues);
        std::cout << "Root of mean squared error: " << rms << ", with a maximum deviation of " << maxDiff << " at ("
                  << maxX << "," << maxY << "," << maxZ << ")[" << fieldMax << "]"
                  << " and average percentile deviation " << averageDeviation * 100.0 << "." << std::endl;
    }
    return averageDeviation;
}

template <unsigned padding>
void GridFunctions::writeToFile(const std::string fileName, const PaddedGrid<FieldStruct, padding>& grid) {
    std::ofstream outputStream;
    outputStream.open(fileName, std::ios::out);
    outputStream.precision(20);
    outputStream << grid.xEnd() - grid.xStart() << " " << grid.yEnd() - grid.yStart() << " "
                 << grid.zEnd() - grid.zStart() << std::endl;
    for (unsigned x = grid.xStart(); x < grid.xEnd(); x++) {
        for (unsigned y = grid.yStart(); y < grid.yEnd(); y++) {
            for (unsigned z = grid.zStart(); z < grid.zEnd(); z++) {
                for (unsigned field = 0; field < NUM_PHYSICAL_FIELDS; field++) {
                    outputStream << std::fixed << grid(x, y, z)[field] << std::endl;
                }
            }
        }
    }
}

namespace GridFunctionsSycl {

bool checkNaN(sycl::queue& queue, sycl::buffer<FieldStruct, 3>& grid);

} // namespace GridFunctionsSycl

namespace GridFunctionsCelerity {

bool checkNaN(celerity::distr_queue& queue, celerity::buffer<FieldStruct, 3>& grid);

} // namespace GridFunctionsCelerity
